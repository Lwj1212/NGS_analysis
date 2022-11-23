# RNA-Seq

---

# RNA-Seq Pipeline

- **Pipeline(Genome mapping or Trascriptome mapping)**

    ![png](https://camo.githubusercontent.com/11feecf7f1ede165c1775cd0ede84a2c9193117df92166eefa991f13f3397dfc/68747470733a2f2f62696f636f72656372672e6769746875622e696f2f524e417365715f636f757273655f323031392f696d616765732f524e417365715f776f726b666c6f772e706e67)

- **Pipeline**
    - Genome mapping - **Tophat2 + Cufflinks, STAR + Cufflinks**
    - Transcriptome mapping - **Salmon, Kallisto, STAR + RSEM**
- **Reference**
    - Genome
        - DNA : [Ensembl GRCh38(hg38)](http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz)
        - GTF : [Ensembl release 104](http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz)
        - cDNA(Transcriptome) : [Ensembl release 104](http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz)
    - **Index**

        [](http://192.168.0.90:5000/sharing/cxH01fTHA)

- **Script**

    ```bash
    git clone https://github.com/Jin0331/RNA_script.git
    ```

    [Jin0331/RNA_script](https://github.com/Jin0331/RNA_script)

    [](https://camo.githubusercontent.com/11feecf7f1ede165c1775cd0ede84a2c9193117df92166eefa991f13f3397dfc/68747470733a2f2f62696f636f72656372672e6769746875622e696f2f524e417365715f636f757273655f323031392f696d616765732f524e417365715f776f726b666c6f772e706e67)

    ## Pre - Processing

    1. **sra_fasatq_downloader.sh**

        ```bash
        # ENA Browser에서 Public data 다운로드
        # SRR 또는 ERR만 가능(single / paired-end 상관없음)
        # ENA Browser에서 받은 파일
        # 다운로드 받을 SRR 또는 ERR의 Accesion Number를 특정 파일에 생성할 것.
        nano sra_input.txt

        RNA_script/script/sra_fastq_downloader.sh **[다운로드받을 폴더의 PATH, ex) /home/sempre813/PUBLIC]** sra_input.txt
        ```

        - sra_input.txt 예시

            [sra_input.txt](RNA-Seq%20abbf1d74591b4d87b227bc929420d5f9/sra_input.txt)

    2. **ena_file_rename.py(ENA Browser에서 받은 fastq 파일(single/paired 모두)만 적용할 것!!)**

        ```bash
        # ENA Browser에서 다운받은 fastq 파일은 single의 경우 _1.fastq.gq, paired의 경우
        # _1.fastq.gz / _2.fastq.gz로 구성 됨.
        # 생성한 Pipeline Script의 경우, _R1.fastq.gz, _R2.fastq.gz로 짜여져 있음.
        # 따라서 _1.fastq.gz / _2.fastqs.gz를 _R1.fastq.gz / _R2.fastq.gz로 변경해야 함.

        python3 RNA_script/script/ena_file_rename.py [*.fastq.gz 파일이 존재하는 디렉토리의 PATH]
        ```

    ## Pipeline 적용

    ENA Browser의 경우 **sra_fasatq_downloader.sh와 ena_file_rename.py**를 이용하여 input 형식에 맞는 파일을 생성할 수 있음. ****하지만~!@ In-House의 경우는 직접 일일히 _R1.fastq.gz, _R2.fastq.gz를 붙여야 함.. 이것도 관련 스크립트를 생성할 계획임.

    ### 1. Genome mapping

    ### 2. Transcriptome mapping

    - **Paired-end**

    **kallisto_salmon.sh**

    ```bash
    ** Argument
        -d : Docker container name
    	-r : a - Salmon, Kallisto / k - Kallisto / s - Salmon
    	**-n : NAS ID - index 폴더가 없는 경우 사용할 것
    	-p : NAS PASSWORD - index 폴더가 없는 경우 사용할 것
    	-h : NAS ip - index 폴더가 없는 경우 사용할 것**
    	-b : base directory - **input과 index 디렉토리가 존재하는 PHAH**
    	-i : input directory - **input 디렉토리의 이름**
    	-t : thread 개수
    	-g : fastq or fq or qc - fq는 *_R1.fq.gz, *_R2.fq.gz 일때만 사용
                             - qc는 trim_galore만 완료되었을 때 사용
        -c : output_dir/QC에 QC된 파일이 있는 경우 y. 
    # 예시
    RNA_script/script/**kallisto_salmon.sh** -r a -b /home/wmbio/ -i PUBLIC_DATA -t 15 -g fastq
    ```

    **star_rsem.sh**

    ```bash
    ** Argument
    	-r : a - STAR+RSEM / s - STAR / r - RSEM
    	**-n : NAS ID - index 폴더가 없는 경우 사용할 것
    	-p : NAS PASSWORD - index 폴더가 없는 경우 사용할 것
    	-h : NAS ip - index 폴더가 없는 경우 사용할 것**
    	-b : base directory - **input과 index 디렉토리가 존재하는 PHAH**
    	-i : input directory - **input 디렉토리의 이름**
    	-t : thread 개수
    	-g : fastq or fq or qc - fq는 *_R1.fq.gz, *_R2.fq.gz 일때만 사용
                             - qc는 trim_galore만 완료되었을 때 사용

    # 예시
    RNA_script/script/**star_rsem.sh** -r a -b /home/wmbio/ -i PUBLIC-2-5 -t 15 -g fastq -s star_ensembl_RON_ADD -e rsem_ensembl_RON_ADD
    ```

- **<Single-end>**

    ```bash
    # 위와 사용법은 동일하며, **kallisto_salmon.sh -> kallisto_salmon_single.sh
    #                   star_rsem.sh -> star_rsem_single.sh
    # 로 변경하여 사용하면 됨.**
    ```
