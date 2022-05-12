version 1.0 

task fastqtl_nominal {
    input {
        File expression_bed
        File expression_bed_index
        File vcf
        File vcf_index
        String prefix
        File covariates

        String? cis_window
        Int? ma_sample_threshold
        Float? maf_threshold
        Int chunks

        Int memory=10
        Int disk_space
        Int num_threads
        Int num_preempt
    }

    command <<<
        set -euo pipefail

        wget https://raw.githubusercontent.com/broadinstitute/palantir-workflows/main/Scripts/monitoring/cromwell_monitoring_script.sh
        bash ./cromwell_monitoring_script.sh | tee monitoring.log &

        touch ~{vcf_index}  # avoid tabix "index older than vcf" error
        touch ~{expression_bed_index}
        # nominal pass

        cat -<< "EOF" > run_FastQTL_threaded.py
        #!/usr/bin/env python3
        # Author: Francois Aguet
        import argparse
        import os
        import numpy as np
        import subprocess
        import gzip
        import multiprocessing as mp
        import contextlib
        from datetime import datetime
        import tempfile
        import glob

        @contextlib.contextmanager
        def cd(cd_path):
            saved_path = os.getcwd()
            os.chdir(cd_path)
            yield
            os.chdir(saved_path)

        def get_cmd(args, chunk):
            cmd = os.path.join(fastqtl_dir, 'bin', 'fastQTL')+' --vcf '+args.vcf+' --bed '+args.bed+' --window '+args.window \
                +' --maf-threshold '+args.maf_threshold \
                +' --ma-sample-threshold '+args.ma_sample_threshold \
                +' --interaction-maf-threshold '+args.interaction_maf_threshold
            if args.covariates:
                cmd += ' --cov '+args.covariates
            if args.phenotype_groups:
                cmd += ' --grp '+args.phenotype_groups
            if args.threshold:
                cmd += ' --threshold '+args.threshold
            if args.permute:
                cmd += ' --permute '+' '.join([str(p) for p in args.permute])
            if args.interaction:
                cmd += ' --interaction '+args.interaction
            if args.best_variant_only:
                cmd += ' --report-best-only'
            if args.seed:
                cmd += ' --seed '+args.seed
            if args.exclude_samples:
                cmd += ' --exclude-samples '+args.exclude_samples
            if args.exclude_sites:
                cmd += ' --exclude-sites '+args.exclude_sites
            cmd += ' --chunk '+str(chunk)+' '+args.chunks\
                + ' --out '+args.prefix+'_chunk{0:03d}.txt.gz'.format(chunk)\
                + ' --log '+args.prefix+'_chunk{0:03d}.log'.format(chunk)
            return cmd

        def perm_worker(inputs):
            args = inputs[0]
            chunk = inputs[1]
            cmd = get_cmd(args, chunk)
            print('Processing chunk '+str(chunk), flush=True)
            s = subprocess.check_call(cmd, shell=True, executable='/bin/bash', stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            print('Finished chunk '+str(chunk), flush=True)
            return s


        parser = argparse.ArgumentParser(description='Run FastQTL')
        parser.add_argument('vcf', help='Genotypes in VCF 4.1 format')
        parser.add_argument('bed', help='Phenotypes in UCSC BED extended format')
        parser.add_argument('prefix', help='Prefix for output file name')
        parser.add_argument('--covariates', default='', help='Covariates')
        parser.add_argument('--phenotype_groups', default='', help='File with mapping of phenotype_id to group_id (gene_id)')
        parser.add_argument('--chunks', default='100', help='Number of chunks, minimum: #chromosomes')
        parser.add_argument('--permute', default=None, type=str, nargs='+', help='Number of permutations, e.g. [1000, 10000] (adaptive). Default: None (run nominal pass)')
        parser.add_argument('--interaction', default=None, type=str, help='Interaction term')
        parser.add_argument('--best_variant_only', action='store_true')
        parser.add_argument('--window', default='1e6', help='Cis-window size. Default values is 1Mb (1e6).')
        parser.add_argument('--threshold', default='', help='Output only significant phenotype-variant pairs with a p-value below threshold (default 1)')
        parser.add_argument('--maf_threshold', default='0.0', help='Include only genotypes with minor allele frequency >=maf_threshold (default 0)')
        parser.add_argument('--ma_sample_threshold', default='0', help='Include only genotypes with >=ma_sample_threshold samples carrying the minor allele (default 0)')
        parser.add_argument('--interaction_maf_threshold', default='0', help='MAF threshold for interactions, applied to lower and upper half of samples')
        parser.add_argument('--fdr', default=0.05, type=np.double)
        parser.add_argument('--seed', default=None, help='Random number generator seed')
        parser.add_argument('--exclude_samples', default=None, help='')
        parser.add_argument('--exclude_sites', default=None, help='')
        parser.add_argument('--qvalue_lambda', default=None, help='lambda parameter for pi0est in qvalue.')
        parser.add_argument('-t', '--threads', default=8, type=int, help='Number of threads')
        parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
        args = parser.parse_args()
        fastqtl_dir =  "/opt/fastqtl"

        if not os.path.exists(args.output_dir):
            os.makedirs(args.output_dir)

        print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Running FastQTL on {0:d} threads.'.format(args.threads), flush=True)

        with cd(args.output_dir):
            with mp.Pool(processes=args.threads) as pool:
                pdata_res = [pool.map_async(perm_worker, ((args,k),)) for k in np.arange(1,int(args.chunks)+1)]
                pool.close()
                pool.join()
            print(pdata_res)

            for res in pdata_res:  # check exit status
                if res is None:
                    print("one of the worker threads returned a None!")
                    assert False
                try:
                    resget=res.get()
                except Exception as e:
                    print("An exception got raised in one of the jobs:")
                    print(res)
                    print(e)
                    raise(e)

                if resget is None:
                    print("one of the worker threads res.get() returned None!")
                    print(res)
                if len(res.get())==0:
                    print("one of the worker threads returned an empty array:")
                    print(res)
                    assert False  
                if res.get()[0]!=0:
                    print("one of the worker threads returned a non-zero value:")
                    print(res)
                    assert False



            with tempfile.NamedTemporaryFile(mode='w+') as chunk_list_file, \
                 tempfile.NamedTemporaryFile(mode='w+') as log_list_file:

                # write chunk and log paths to file
                chunk_files = sorted(glob.glob(args.prefix+'_chunk*.txt.gz'))
                chunk_list_file.write('\n'.join(chunk_files)+'\n')
                chunk_list_file.flush()
                log_files = sorted(glob.glob(args.prefix+'_chunk*.log'))
                log_list_file.write('\n'.join(log_files)+'\n')
                log_list_file.flush()

                # merge chunks
                cmd = 'python3 '+os.path.join(fastqtl_dir, 'python', 'merge_chunks.py') \
                    +' {} {} {} --fdr {} -o .'.format(chunk_list_file.name, log_list_file.name, args.prefix, args.fdr)
                if args.qvalue_lambda:
                    cmd += ' --qvalue_lambda {}'.format(args.qvalue_lambda)
                if args.permute:
                    cmd += ' --permute'
                subprocess.check_call(cmd, shell=True)

                # remove chunk files
                for f in chunk_files + log_files:
                    os.remove(f)
        EOF

        env python3 run_FastQTL_threaded.py ~{vcf} ~{expression_bed} ~{prefix} \
            --covariates ~{covariates} \
            ~{"--window " + cis_window} \
            ~{"--ma_sample_threshold " + ma_sample_threshold} \
            ~{"--maf_threshold " + maf_threshold} \
            --chunks ~{chunks} \
            --threads ~{num_threads}

    >>>

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_threads}"
        preemptible: "~{num_preempt}"
    }

    output {
        Array[File] logs = glob("*.log")
        File monitoring_log="monitoring.log"
        File allpairs_log="~{prefix}.allpairs.log"
        File? allpairs="~{prefix}.allpairs.txt.gz" #optional so that delocalization doesn't stop if missing
    }

    meta {
        author: "Francois Aguet"
    }
}


task fastqtl_permutations_scatter {
    input {
        File expression_bed
        File expression_bed_index
        File vcf
        File vcf_index
        String prefix
        File covariates

        Int current_chunk
        Int chunks
        Int permutations
        String? cis_window
        File? phenotype_groups
        Int? ma_sample_threshold
        Float? maf_threshold

        Int memory
        Int disk_space
        Int num_threads
        Int num_preempt
    }

    command <<<
        set -euo pipefail
        touch ~{vcf_index}  # avoid tabix "index older than vcf" error
        touch ~{expression_bed_index}
        # permutation pass
        /opt/fastqtl/python/run_chunk.py ~{vcf} ~{expression_bed} ~{prefix} ~{current_chunk} ~{chunks}\
            --permute ~{permutations} \
            --covariates ~{covariates} \
            ~{"--window " + cis_window} \
            ~{"--phenotype_groups " + phenotype_groups} \
            ~{"--ma_sample_threshold " + ma_sample_threshold} \
            ~{"--maf_threshold " + maf_threshold}
        mv ~{prefix}_chunk*.txt.gz ~{prefix}_chunk_~{current_chunk}.txt.gz
        mv ~{prefix}_chunk*.log ~{prefix}_chunk_~{current_chunk}.log
    >>>

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_threads}"
        preemptible: "~{num_preempt}"
    }

    output {
        File chunk="~{prefix}_chunk_~{current_chunk}.txt.gz"
        File chunk_log="~{prefix}_chunk_~{current_chunk}.log"
    }

    meta {
        author: "Francois Aguet"
    }
}


task fastqtl_permutations_merge {
    input {
        Array[File] chunks
        Array[File] logs
        String prefix

        Int memory
        Int disk_space
        Int num_threads
        Int num_preempt

        Float? qvalue_lambda
    }
    command <<<
        set -euo pipefail
        /opt/fastqtl/python/merge_chunks.py ~{write_lines(chunks)} ~{write_lines(logs)} ~{prefix}\
            --permute ~{"--qvalue_lambda " + qvalue_lambda}
    >>>

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_threads}"
        preemptible: "~{num_preempt}"
    }

    output {
        File genes="~{prefix}.genes.txt.gz"
        File genes_log="~{prefix}.genes.log"
    }

    meta {
        author: "Francois Aguet"
    }
}


task fastqtl_postprocess {
    input {
        File permutations_output
        File nominal_output
        Float fdr
        File annotation_gtf
        String prefix
        File? variant_lookup

        Int memory
        Int disk_space
        Int num_threads
        Int num_preempt
    }
    command <<<
        set -euo pipefail

        cat -<< "EOF" > annotate_outputs.py
        #!/usr/bin/env python3
        # Author: Francois Aguet

        import argparse
        import numpy as np
        import pandas as pd
        import os
        import gzip
        from datetime import datetime
        import subprocess
        import io

        parser = argparse.ArgumentParser(description='Filter significant SNP-gene pairs from FastQTL results using FDR cutoff')
        parser.add_argument('permutation_results', help='FastQTL output')
        parser.add_argument('fdr', type=np.double, help='False discovery rate (e.g., 0.05)')
        parser.add_argument('annotation_gtf', help='Annotation in GTF format')
        parser.add_argument('--snp_lookup', default='', help='Tab-delimited file with columns: chr, variant_pos, variant_id, ref, alt, num_alt_per_site, rs_id_dbSNP...')
        parser.add_argument('--nominal_results', default='', help='FastQTL output from nominal pass')
        parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
        args = parser.parse_args()

        #------------------------------------------------------------------------------
        # 1. eGenes (permutation output): add gene and variant information
        #------------------------------------------------------------------------------
        gene_dict = {}
        print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Parsing GTF', flush=True)
        # add: gene_name, gene_chr, gene_start, gene_end, strand
        with gzip.open(args.annotation_gtf,'rt') as gtf:
            for row in gtf:
                row = row.strip().split('\t')
                if row[0][0]=='#' or row[2]!='gene': continue
                # get gene_id and gene_name from attributes
                attr = dict([i.split() for i in row[8].replace('"','').split(';') if i!=''])
                gene_dict[attr['gene_id']] = [attr['gene_name'], row[0], row[3], row[4], row[6]]

        print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Annotating permutation results (eGenes)', flush=True)
        gene_df = pd.read_csv(args.permutation_results, sep='\t', index_col=0)
        if 'group_id' in gene_df:
            gene_info = pd.DataFrame(data=[gene_dict[i] for i in gene_df['group_id']], columns=['gene_name', 'gene_chr', 'gene_start', 'gene_end', 'strand'], index=gene_df.index)
        else:
            attr = {i:i.split(":")[4] for i in gene_df.index}
            gene_info = pd.DataFrame(data=[gene_dict[attr[i]] for i in gene_df.index], columns=['gene_name', 'gene_chr', 'gene_start', 'gene_end', 'strand'], index=gene_df.index)
            
        gene_df = pd.concat([gene_info, gene_df], axis=1)
        assert np.all(gene_df.index==gene_info.index)

        col_order = ['gene_name', 'gene_chr', 'gene_start', 'gene_end', 'strand',
            'num_var', 'beta_shape1', 'beta_shape2', 'true_df', 'pval_true_df', 'variant_id', 'tss_distance']
        if args.snp_lookup:
            print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Adding variant annotations from lookup table', flush=True)
            # intersect lookup table with variant_ids (col 7 in permutation_results; col 3 in snp_lookup)
            cmd = "awk 'NR==FNR {v[$7]; next} $3 in v' <(zcat " + args.permutation_results + ") <(zcat " + args.snp_lookup + ")"
            s = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
            snp_lookup_df = pd.read_csv(io.StringIO(s.decode()), index_col=2, sep='\t',
                dtype={'chr':str, 'variant_pos':np.int64, 'variant_id':str, 'ref':str, 'alt':str, 'num_alt_per_site':np.int32})
            gene_df = gene_df.join(snp_lookup_df, on='variant_id')  # add variant information
            col_order += list(snp_lookup_df.columns)
        col_order += ['ma_samples', 'ma_count', 'maf', 'ref_factor',
            'pval_nominal', 'slope', 'slope_se', 'pval_perm', 'pval_beta']
        if 'group_id' in gene_df:
            col_order += ['group_id', 'group_size']
        col_order += ['qval', 'pval_nominal_threshold']
        gene_df = gene_df[col_order]

        outname = os.path.join(args.output_dir, os.path.split(args.permutation_results)[1].split('.txt.gz')[0]+'.annotated.txt.gz')
        with gzip.open(outname, 'wt') as f:
            gene_df.to_csv(f, sep='\t', float_format='%.6g')

        #------------------------------------------------------------------------------
        # 2. variant-gene pairs: output new file with all significant pairs
        #------------------------------------------------------------------------------
        if args.nominal_results:
            print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Filtering significant variant-gene pairs', flush=True)

            # eGenes (apply FDR threshold)
            egene_df = gene_df.loc[gene_df['qval']<=args.fdr, ['pval_nominal_threshold', 'pval_nominal', 'pval_beta']].copy()
            egene_df.rename(columns={'pval_nominal': 'min_pval_nominal'}, inplace=True)
            egene_ids = set(egene_df.index)
            threshold_dict = egene_df['pval_nominal_threshold'].to_dict()

            # process by chunks to reduce memory usage
            signif_df = []
            mask = []
            for i,chunk in enumerate(pd.read_csv(args.nominal_results, sep='\t', iterator=True, chunksize=1000000, index_col=1,
                dtype={'gene_id':str, 'variant_id':str, 'tss_distance':np.int32,
                    'ma_samples':np.int32, 'ma_count':np.int32, 'maf':np.float32,
                    'pval_nominal':np.float64, 'slope':np.float32, 'slope_se':np.float32})):
                chunk = chunk[chunk['gene_id'].isin(egene_ids)]
                m = chunk['pval_nominal']<chunk['gene_id'].apply(lambda x: threshold_dict[x])
                signif_df.append(chunk[m])
                mask.append(m)
                print('Chunks processed: {0:d}'.format(i+1), end='\r', flush=True)
            signif_df = pd.concat(signif_df, axis=0)
            signif_df = signif_df.merge(egene_df, left_on='gene_id', right_index=True)

            outname = os.path.join(args.output_dir, os.path.split(args.nominal_results)[1].split('.allpairs.txt.gz')[0]+'.signifpairs.txt.gz')
            with gzip.open(outname, 'wt') as f:
                signif_df.to_csv(f, sep='\t', float_format='%.6g')

        print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Completed annotation', flush=True)

        EOF

        # post-processing
        python3 ./annotate_outputs.py \
            ~{permutations_output} \
            ~{fdr} \
            ~{annotation_gtf} \
            --nominal_results ~{nominal_output} \
            ~{"--snp_lookup " + variant_lookup}
    >>>

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_threads}"
        preemptible: "~{num_preempt}"
    }

    output {
        File genes_annotated="~{prefix}.genes.annotated.txt.gz"
        File signifpairs="~{prefix}.signifpairs.txt.gz"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow fastqtl_workflow {
    input {
        File expression_bed
        File expression_bed_index
        File vcf
        File vcf_index
        String prefix
        File covariates

        Int permutations
        Int chunks
        String? cis_window
        Int? ma_sample_threshold
        Float? maf_threshold

        # post-processing
        Float fdr
        File annotation_gtf
        File? variant_lookup
    }
    call fastqtl_nominal {
        input:
            chunks=chunks, 
            prefix=prefix,
            expression_bed=expression_bed, 
            expression_bed_index=expression_bed_index, 
            vcf=vcf, 
            vcf_index=vcf_index, 
            covariates=covariates, 
            cis_window=cis_window,
            ma_sample_threshold=ma_sample_threshold, 
            maf_threshold=maf_threshold,
            disk_space=ceil(3*size(vcf,"GB")+200),
            num_threads=4,
            num_preempt=1
    }

    scatter(i in range(chunks)) {
        call fastqtl_permutations_scatter {
            input:
                current_chunk=i+1, 
                chunks=chunks, 
                prefix=prefix, 
                permutations=permutations,
                expression_bed=expression_bed, 
                expression_bed_index=expression_bed_index,
                vcf=vcf, 
                vcf_index=vcf_index,
                covariates=covariates, 
                cis_window=cis_window,
                ma_sample_threshold=ma_sample_threshold, 
                maf_threshold=maf_threshold,
                memory=10,
                disk_space=ceil(size(vcf,"GB")+20),
                num_threads=4,
                num_preempt=1
        }
    }

    call fastqtl_permutations_merge {
        input: 
            chunks=fastqtl_permutations_scatter.chunk, 
            logs=fastqtl_permutations_scatter.chunk_log, 
            prefix=prefix,
            memory=10,
            disk_space=50,
            num_threads=4,
            num_preempt=1
        }

    output {
        File genes_permutation=fastqtl_permutations_merge.genes
        File? genes_nominal=fastqtl_nominal.allpairs

   }
}
