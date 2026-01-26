nextflow.enable.dsl=2

/*
 * PARAMETERS
 */
params.input         = 'data'
params.outdir        = 'results'
params.ligand_results = 'results/ligands'
params.prep_receptor = "${projectDir}/envs/slc6env/bin/prepare_receptor4.py"




process RUN_FPOCKET {

    tag "fpocket"

    publishDir params.outdir, mode: 'copy'

    cpus 2
    memory '8 GB'
    time '2h'

    input:
    path input_dir

    output:
    path "fpocket_results"

    script:
    """
    mkdir fpocket_results

    ${projectDir}/src/fpocket.sh ${input_dir}

    # fpocket crée *_out → on les déplace dans le dossier attendu
    mv *_out fpocket_results/ 2>/dev/null || true
    """
}


/*
 * PROCESS 2: POCKET SCORING
 */
process POCKET_SCORING {

    tag "pocket_scoring"


    publishDir params.outdir, mode: 'copy'

    cpus 2
    memory '4 GB'
    time '1h'

    input:
    path template_pdb
    path target_pdb
    path fpocket_results

    output:
    path "pocket_scores_ranked.csv"

    script:
    """
    python3 ${projectDir}/src/pocket_scoring.py \
        ${template_pdb} \
        ${target_pdb} \
        ${fpocket_results}
    """
}

/*
 * PROCESS 3: RESOLVE SEEDS (ChEMBL / SMILES)
 */
process RESOLVE_SEEDS {

    tag "resolve_seeds"

    publishDir "${params.ligand_results}/seeds", mode: 'copy'

    input:
    path seeds_tsv

    output:
    path "seeds.sdf"
    path "seeds_resolve.tsv"

    script:
    """
    python3 ${projectDir}/src/resolve_seeds_chembl.py \
        --in_tsv ${seeds_tsv} \
        --out_sdf seeds.sdf \
        --out_report_tsv seeds_resolve.tsv
    """
}


/*
 * PROCESS 4: ChEMBL HOMOLOG SEARCH
 */
process CHEMBL_HOMOLOGS {

    tag "chembl_homologs"

    publishDir "${params.ligand_results}/homologs", mode: 'copy'

    cpus 1
    memory '4 GB'
    time '2h'

    input:
    path seeds_sdf

    output:
    path "chembl_homologs.sdf"
    path "chembl_edges.tsv"

    script:
    """
    python3 ${projectDir}/src/chembl_homolog_search.py \
        --in_sdf ${seeds_sdf} \
        --out_sdf chembl_homologs.sdf \
        --out_edges chembl_edges.tsv
    """
}


/*
 * PROCESS 5: LIGANDS TO 3D
 */
process LIGANDS_3D {

    tag "ligands_3d"

    publishDir "${params.ligand_results}/3d", mode: 'copy'

    input:
    path homologs_sdf

    output:
    path "ligands_3d.sdf"

    script:
    """
    python3 ${projectDir}/src/ligands_to_3d.py \
        --in_sdf ${homologs_sdf} \
        --out_sdf ligands_3d.sdf
    """
}

/*
 * PROCESS 6: MINIMIZE LIGANDS
 */
process MINIMIZE_LIGANDS {

    tag "minimize"

    publishDir "${params.ligand_results}/minimized", mode: 'copy'

    input:
    path ligands_3d_sdf

    output:
    path "ligands_minimized.sdf"

    script:
    """
    python3 ${projectDir}/src/minimize_ligands.py \
        --in_sdf ${ligands_3d_sdf} \
        --out_sdf ligands_minimized.sdf
    """
}

/*
 * PROCESS 7: SDF to PDBQT
 */
process SDF_TO_PDBQT {

    tag "pdbqt"

    publishDir "${params.ligand_results}", mode: 'copy'

    input:
    path ligands_min_sdf

    output:
    path "pdbqt"


    script:
    """
    python3 ${projectDir}/src/sdf_to_pdbqt_obabel.py \
        --in_sdf ${ligands_min_sdf} \
        --out_dir pdbqt
    """
}

/*
 * PROCESS 8: DOCKING
 */
process DOCKING {
    tag "vina_docking"
    publishDir params.outdir, mode: 'copy'
    cpus params.docking_cpu
    memory '8 GB'
    time '12h'
    beforeScript 'export PATH="/shared/ifbstor1/projects/tp_2549_meet_eu_2526_182839/su_team_2/envs/slc6env/bin:$PATH"'
    input:
    path template_pdb
    path target_pdb
    path ligands_pdbqt_dir
    path pocket_scores_csv
    output:
    path "docking_results"
    script:
    """
    python3 ${projectDir}/src/docking.py \
        --template_pdb ${template_pdb} \
        --target_pdb ${target_pdb} \
        --ligands ${ligands_pdbqt_dir} \
        --pocket_scores_csv ${pocket_scores_csv} \
        --out docking_results \
        --prep_receptor ${params.prep_receptor}
    """
}

/*
 * WORKFLOW
 */
workflow {

    // Structure / Pocket part
    input_dir    = file(params.input)
    template_pdb = file("data/template/positive/8WBY_chainA.pdb")
    target_pdb   = file("data/target/positive/SIT1_model_oo.pdb")

    fpocket_out = RUN_FPOCKET(input_dir)
    pocket_csv  = POCKET_SCORING(template_pdb, target_pdb, fpocket_out)

    // Ligand / ChEMBL part
    seeds_tsv = file("data/ligands/seeds/known_inhibitors.tsv")

    seeds    = RESOLVE_SEEDS(seeds_tsv)
    homologs = CHEMBL_HOMOLOGS(seeds[0])
    lig3d    = LIGANDS_3D(homologs[0])
    ligmin   = MINIMIZE_LIGANDS(lig3d)

    ligands_pdbqt_dir = SDF_TO_PDBQT(ligmin)
    

    // Docking step
    DOCKING(
        template_pdb,
        target_pdb,
        ligands_pdbqt_dir,
        pocket_csv
    )
}
