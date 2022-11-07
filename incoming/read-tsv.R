library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
read.files <- function()
{
    files <- list.files(pattern=".txt")

    # [1] "CE-only.elib.proteins.txt"
    # [2] "NE1_Prosit_based_quant report.elib.proteins.txt"
    # [3] "NE2_DIA_quant_filtered_10252022.txt"



    #------------------------------------------------
    # cyto: read and remove multiply-mapped proteins
    #------------------------------------------------

    tbl.cyto <-
        read.table(files[1], sep="\t", header=TRUE, nrow=-1)[,-3]
    dim(tbl.cyto)  # 4806   11
    plural.mapped <- grep(";", tbl.cyto$Protein)
    length(plural.mapped)  # 75
    tbl.cyto <- tbl.cyto[-plural.mapped,]
    dim(tbl.cyto)  # 4731 11

    # identify and remove spiked-in proteins
    non.uniref <- which(!grepl("UniRef90", tbl.cyto$Protein))
    tbl.cyto$Protein[non.uniref]

     # [1] "ZcRAP|ANXA5_HUMAN|" "ZcRAP|TRFE_HUMAN|"  "ZcRAP|ALBU_BOVIN|"  "ZcRAP|CASK_BOVIN|"
     # [5] "ZcRAP|GFP_AEQVI|"   "ZcRAP|CAS1_BOVIN|"  "ZcRAP|TRYP_PIG|"    "ZcRAP|CAS2_BOVIN|"
     # [9] "ZcRAP|REF_HEVBR|"

    tbl.cyto <- tbl.cyto[-non.uniref,]
    dim(tbl.cyto)  #  4722   11
    tbl.cyto$Protein <- sub("UniRef90_", "", tbl.cyto$Protein)
    head(tbl.cyto)

    #   Protein NumPeptides        D0         D2        D4         D6        D8       D10        D11        D12       D14
    # 1  P04264          26 489000000 1620000000 598000000 4800000000 968000000 884000000 1340000000 3480000000 494000000
    # 2  Q92609          13  87700000   49700000  45000000   43600000  46800000  55100000   42500000   39700000  81900000
    # 3  Q99571           3   3830394    4336378   3026220    5367845   5742663   8499565    7710952    6140192   3398962
    # 4  Q8IV04           2  11100000    4521368   5673966    4534420   5836556   4804094    5388049    1307840  10700000
    # 5  P62136           4 129000000  110000000 175000000  158000000 159000000 162000000  157000000  134000000 124000000
    # 6  Q5JSH3          13  44200000   25500000  44400000   49800000  57600000  62800000   63700000   40800000  45200000

    #-----------------------------------------------
    # ne1: read and remove multiply-mapped proteins
    #-----------------------------------------------

    tbl.ne1 <-
        read.table(files[2], sep="\t", header=TRUE, nrow=-1)[, -3]
    dim(tbl.ne1)  # 4049 11
    colnames(tbl.ne1)
    plural.mapped <- grep(";", tbl.ne1$Protein)
    length(plural.mapped)  # 84
    tbl.ne1 <- tbl.ne1[-plural.mapped,]
    # identify and remove spiked-in proteins

    non.uniref <- which(!grepl("UniRef90", tbl.ne1$Protein))
    tbl.ne1$Protein[non.uniref]

      # [1] "ZcRAP|GFP_AEQVI|"   "ZcRAP|ALBU_BOVIN|"  "ZcRAP|ANXA5_HUMAN|"
      # [4] "ZcRAP|TRFE_HUMAN|"  "ZcRAP|REF_HEVBR|"   "ZcRAP|CAS1_BOVIN|"
      # [7] "ZcRAP|TRYP_PIG|"

    tbl.ne1 <- tbl.ne1[-non.uniref,]
    dim(tbl.ne1)  # 3958   11
    tbl.ne1$Protein <- sub("UniRef90_", "", tbl.ne1$Protein)
    head(tbl.ne1)

       # all data columns have same verbose prefix. remove them)
    ncol(tbl.ne1) # 11
    prefix <- "X2_20220822_LM_Eclipse.ES.15cm.300nlml_90min_OTOT_DIA_all_NE1_"
    length(grep(prefix, colnames(tbl.ne1)))
    new.col.names <- sub(prefix, "", colnames(tbl.ne1))
    suffix <- ".mzML"
    new.col.names <- sub(suffix, "", new.col.names, fixed=TRUE)
    new.col.names.timeOrdered <-
        c("Protein", "NumPeptides", "D0", "D2", "D4", "D6", "D8", "D10", "D11", "D12", "D14")

    colnames(tbl.ne1) <- new.col.names.timeOrdered

    #-----------------------------------------------
    # ne2: read and remove multiply-mapped proteins
    #-----------------------------------------------

    tbl.ne2 <-
        read.table(files[3], sep="\t", header=TRUE, nrow=-1)[, -3]
    dim(tbl.ne2)  #  3067 11
    plural.mapped <- grep(";", tbl.ne2$Protein)
    length(plural.mapped)  # 28
    tbl.ne2 <- tbl.ne2[-plural.mapped,]
    dim(tbl.ne2)  # 3000 11
       # identify and remove spiked-in proteins
    non.uniref <- which(!grepl("UniRef90", tbl.ne2$Protein))
    tbl.ne2$Protein[non.uniref]
    tbl.ne2 <- tbl.ne2[-non.uniref,]
    dim(tbl.ne2)  # 2990   11
    tbl.ne2$Protein <- sub("UniRef90_", "", tbl.ne2$Protein)
    head(tbl.ne2)

       #  shorten the colnames
    prefix <- "X3_20220822_LM_Eclipse.ES.15cm.300nlml_90min_OTOT_DIA_all_NE2_"
    length(grep(prefix, colnames(tbl.ne2))) # 9

    new.col.names <- sub(prefix, "", colnames(tbl.ne2))
    suffix <- ".mzML"
    new.col.names <- sub(suffix, "", new.col.names, fixed=TRUE)
    new.col.names.timeOrdered <-
        c("Protein", "NumPeptides", "D0", "D2", "D4", "D6", "D8", "D10", "D11", "D12", "D14")

    colnames(tbl.ne2) <- new.col.names.timeOrdered

    #---------------------------------------------
    # summarize
    #---------------------------------------------

    dim(tbl.cyto) #  4722   11
    dim(tbl.ne1)  #  3958   11
    dim(tbl.ne2)  #  2990   11

    save(tbl.cyto, tbl.ne1, tbl.ne2, file="tbls.after.cleanup.RData")

} # read.files
#----------------------------------------------------------------------------------------------------
lookup.uniprotids.biomart <- function()
{
    require(biomaRt)

    print(load("tbls.after.cleanup.RData"))

    mart <- useMart("ensembl","hsapiens_gene_ensembl")
    tbl.mart <- getBM(c("hgnc_symbol", "gene_biotype", "uniprotswissprot"), mart = mart)
    ids <- unique(c(tbl.ne1$Protein, tbl.ne2$Protein, tbl.cyto$Protein))
    head(tbl.mart)

      #   hgnc_symbol   gene_biotype uniprotswissprot
      # 1      MT-ND1 protein_coding           P03886
      # 2      MT-ND2 protein_coding           P03891
      # 3      MT-CO1 protein_coding           P00395
      # 4      MT-CO2 protein_coding           P00403
      # 5     MT-ATP8 protein_coding           P03928
      # 6     MT-ATP6 protein_coding           P00846

    length(ids) # 5790
    length(intersect(ids, tbl.mart$uniprotswissprot))  # 5670
    missing.uniprots <- setdiff(ids, tbl.mart$uniprotswissprot)
    length(missing.uniprots) # 120

    save(ids, missing.uniprots, tbl.mart, file="tbl.uniprots.mapped.to.geneSymbol.RData")

} # lookup.uniprotids.biomart
#----------------------------------------------------------------------------------------------------
lookup.unmapped.uniprotids <- function()
{
  require(rentrez)
  print(load("tbl.uniprots.mapped.to.geneSymbol.RData"))

  # entrez_fetch(db="protein", id="P63005", rettype="text")
  #  [1] "LOCUS       LIS1_MOUSE               410 aa            linear   ROD 12-OCT-2022\nDEFINITION  RecName: Full=Platelet-activating factor acetylhydrolase IB subunit\n            beta; AltName: Full=Lissencephaly-1 protein; Short=LIS-1; AltName:\n            Full=PAF acetylhydrolase 45 kDa subunit; Short=PAF-AH 45 kDa\n            subunit; AltName: Full=PAF-AH alpha; Short=PAFAH alpha.\nACCESSION   P63005\nVERSION     P63005.2\nDBSOURCE    UniProtKB: locus LIS1_MOUSE, accession P63005;\n            class: standard.\n            extra accessions:O35592,P43035,P81692,Q9R2A6\n            created: Aug 31, 2004.\n            sequence updated: Jan 23, 2007.\n            annotation updated: Oct 12, 2022.\n            xrefs: U95116.1, AAC04610.1, U95120.1, AAC63099.1, U95119.1,\n            AAC63098.1, L25109.1, AAD23059.1, AY189217.1, AAO41716.1,\n            AY189218.1, AAO41717.1, BC014831.1, AAH14831.1, BC026141.1,\n            AAH26141.1, NP_038653.1, 1UUJ_A, 1UUJ_B, 1UUJ_C, 1UUJ_D, 1VYH_C,\n            1VYH_D, 1VYH_G, 1VYH_H, 1VYH_K, 1VYH_L, 1VYH_O, 1VYH_P, 1VYH_S,\n            1VYH_T\n            xrefs (non-sequence databases): CCDS:CCDS25035.1, PDBsum:1UUJ,\n            PDBsum:1VYH, AlphaFoldDB:P63005, SMR:P63005, BioGRID:202013,\n            DIP:DIP-29555N, IntAct:P63005, MINT:P63005,\n            STRING:10090.ENSMUSP00000021091, iPTMnet:P63005,\n            PhosphoSitePlus:P63005, SwissPalm:P63005,\n            REPRODUCTION-2DPAGE:P63005, EPD:P63005, jPOST:P63005, MaxQB:P63005,\n            PaxDb:P63005, PeptideAtlas:P63005, PRIDE:P63005,\n            ProteomicsDB:292334, ProteomicsDB:292335, Antibodypedia:3850,\n            DNASU:18472, Ensembl:ENSMUST00000021091,\n            Ensembl:ENSMUSP00000021091, Ensembl:ENSMUSG00000020745,\n            Ensembl:ENSMUST00000102520, Ensembl:ENSMUSP00000099578,\n            GeneID:18472, KEGG:mmu:18472, UCSC:uc007kce.2, CTD:5048,\n            MGI:109520, VEuPathDB:HostDB:ENSMUSG00000020745, eggNOG:KOG0295,\n            GeneTree:ENSGT00940000155039, HOGENOM:CLU_000288_57_15_1,\n            InParanoid:P63005, OMA:TQECKCV, PhylomeDB:P63005, TreeFam:TF105741,\n            Reactome:R-MMU-141444, Reactome:R-MMU-2467813,\n            Reactome:R-MMU-2500257, Reactome:R-MMU-2565942,\n            Reactome:R-MMU-380259, Reactome:R-MMU-380270,\n            Reactome:R-MMU-380284, Reactome:R-MMU-380320,\n            Reactome:R-MMU-5620912, Reactome:R-MMU-5663220,\n            Reactome:R-MMU-6811436, Reactome:R-MMU-68877,\n            Reactome:R-MMU-8854518, Reactome:R-MMU-9648025, BioGRID-ORCS:18472,\n            ChiTaRS:Pafah1b1, EvolutionaryTrace:P63005, PRO:PR:P63005,\n            Proteomes:UP000000589, RNAct:P63005, Bgee:ENSMUSG00000020745,\n            ExpressionAtlas:P63005, Genevisible:P63005, GO:0008247, GO:0000235,\n            GO:0030424, GO:1904115, GO:0005938, GO:0031252, GO:0090724,\n            GO:0005813, GO:0005737, GO:0005881, GO:0005829, GO:0030426,\n            GO:0005871, GO:0000776, GO:0005875, GO:0015630, GO:0031514,\n            GO:0043005, GO:0043025, GO:0005635, GO:0031965, GO:0048471,\n            GO:0032420, GO:0045202, GO:0031982, GO:0070840, GO:0045505,\n            GO:0042802, GO:0008017, GO:0051010, GO:0051219, GO:0046982,\n            GO:0044877, GO:0001675, GO:0030036, GO:0008344, GO:0001667,\n            GO:0060117, GO:0048854, GO:0016477, GO:0021987, GO:0021895,\n            GO:0007268, GO:0090102, GO:0021540, GO:0043622, GO:0051660,\n            GO:0051649, GO:0000132, GO:0042249, GO:0007281, GO:0021766,\n            GO:1904936, GO:0007254, GO:0021819, GO:0007611, GO:0016042,\n            GO:0051661, GO:0000226, GO:0090176, GO:0031023, GO:0051012,\n            GO:0007017, GO:0097529, GO:0046329, GO:0010977, GO:0007405,\n            GO:0050885, GO:0001764, GO:0051081, GO:0007097, GO:0036035,\n            GO:0045773, GO:0051130, GO:0001961, GO:0061003, GO:0040019,\n            GO:0045931, GO:0009306, GO:0140650, GO:0038026, GO:0043087,\n            GO:0070507, GO:0008090, GO:0017145, GO:0019226, GO:0047496,\n            Gene3D:2.130.10.10, HAMAP:MF_03141, InterPro:IPR017252,\n            InterPro:IPR020472, InterPro:IPR037190, InterPro:IPR006594,\n            InterPro:IPR015943, InterPro:IPR001680, InterPro:IPR019775,\n            InterPro:IPR036322, Pfam:PF08513, Pfam:PF00400, PIRSF:PIRSF037647,\n            PRINTS:PR00320, SMART:SM00667, SMART:SM00320, SUPFAM:SSF109925,\n            SUPFAM:SSF50978, PROSITE:PS50896, PROSITE:PS00678, PROSITE:PS50082,\n            PROSITE:PS50294\nKEYWORDS    3D-structure; Acetylation; Alternative splicing; Cell cycle; Cell\n            division; Coiled coil; Cytoplasm; Cytoskeleton; Developmental\n            protein; Differentiation; Direct protein sequencing; Lipid\n            degradation; Lipid metabolism; Membrane; Microtubule; Mitosis;\n            Neurogenesis; Nucleus; Phosphoprotein; Reference proteome; Repeat;\n            Transport; WD repeat.\nSOURCE      Mus musculus (house mouse)\n  ORGANISM  Mus musculus\n            Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;\n            Mammalia; Eutheria; Euarchontoglires; Glires; Rodentia; Myomorpha;\n            Muroidea; Muridae; Murinae; Mus; Mus.\nREFERENCE   1  (residues 1 to 410)\n  AUTHORS   Peterfy,M., Gyuris,T., Basu,R. and Takacs,L.\n  TITLE     Lissencephaly-1 is one of the most conserved proteins between mouse\n            and human: a single amino-acid difference in 410 residues\n  JOURNAL   Gene 150 (2), 415-416 (1994)\n   PUBMED   7821822\n  REMARK    NUCLEOTIDE SEQUENCE [MRNA] (ISOFORM 1).;\n            STRAIN=ICR; TISSUE=Brain\nREFERENCE   2  (residues 1 to 410)\n  AUTHORS   Peterfy,M., Gyuris,T., Grosshans,D., Cuaresma,C.C. and Takacs,L.\n  TITLE     Cloning and characterization of cDNAs and the gene encoding the\n            mouse platelet-activating factor acetylhydrolase Ib alpha\n            subunit/lissencephaly-1 protein\n  JOURNAL   Genomics 47 (2), 200-206 (1998)\n   PUBMED   9479492\n  REMARK    NUCLEOTIDE SEQUENCE [GENOMIC DNA], AND ALTERNATIVE SPLICING.;\n            TISSUE=Brain\nREFERENCE   3  (residues 1 to 410)\n  AUTHORS   Reiner,O., Albrecht,U., Gordon,M., Chianese,K.A., Carrozzo,R.,\n            Wong,C., Lindsay,E.A., Siracusa,L.D., Baldini,A., Ledbetter,D.H.,\n            Eichele,G., Buchberg,A.M. and Caskey,T.C.\n  TITLE     Lissencephaly gene (LIS1) expression in the CNS suggests a role in\n            neuronal migration\n  JOURNAL   J Neurosci 15 (5 Pt 2), 3730-3738 (1995)\n   PUBMED   7751941\n  REMARK    NUCLEOTIDE SEQUENCE [MRNA] OF 115-410 (ISOFORM 2).;\n            TISSUE=Brain\nREFERENCE   4  (residues 1 to 410)\n  AUTHORS   Grattan,M., Mi,Q.S., Meagher,C. and Delovitch,T.L.\n  TITLE     Congenic mapping of the diabetogenic locus Idd4 to a 5.2-cM region\n            of chromosome 11 in NOD mice: identification of two potential\n            candidate subloci\n  JOURNAL   Diabetes 51 (1), 215-223 (2002)\n   PUBMED   11756344\n  REMARK    NUCLEOTIDE SEQUENCE [MRNA] (ISOFORM 1).;\n            STRAIN=C57BL/6J, and NOD\nREFERENCE   5  (residues 1 to 410)\n  AUTHORS   Gerhard,D.S., Wagner,L., Feingold,E.A., Shenmen,C.M., Grouse,L.H.,\n            Schuler,G., Klein,S.L., Old,S., Rasooly,R., Good,P., Guyer,M.,\n            Peck,A.M., Derge,J.G., Lipman,D., Collins,F.S., Jang,W., Sherry,S.,\n            Feolo,M., Misquitta,L., Lee,E., Rotmistrovsky,K., Greenhut,S.F.,\n            Schaefer,C.F., Buetow,K., Bonner,T.I., Haussler,D., Kent,J.,\n            Kiekhaus,M., Furey,T., Brent,M., Prange,C., Schreiber,K.,\n            Shapiro,N., Bhat,N.K., Hopkins,R.F., Hsie,F., Driscoll,T.,\n            Soares,M.B., Casavant,T.L., Scheetz,T.E., Brown-stein,M.J.,\n            Usdin,T.B., Toshiyuki,S., Carninci,P., Piao,Y., Dudekula,D.B.,\n            Ko,M.S., Kawakami,K., Suzuki,Y., Sugano,S., Gruber,C.E.,\n            Smith,M.R., Simmons,B., Moore,T., Waterman,R., Johnson,S.L.,\n            Ruan,Y., Wei,C.L., Mathavan,S., Gunaratne,P.H., Wu,J., Garcia,A.M.,\n            Hulyk,S.W., Fuh,E., Yuan,Y., Sneed,A., Kowis,C., Hodgson,A.,\n            Muzny,D.M., McPherson,J., Gibbs,R.A., Fahey,J., Helton,E.,\n            Ketteman,M., Madan,A., Rodrigues,S., Sanchez,A., Whiting,M.,\n            Madari,A., Young,A.C., Wetherby,K.D., Granite,S.J., Kwong,P.N.,\n            Brinkley,C.P., Pearson,R.L., Bouffard,G.G., Blakesly,R.W.,\n            Green,E.D., Dickson,M.C., Rodriguez,A.C., Grimwood,J., Schmutz,J.,\n            Myers,R.M., Butterfield,Y.S., Griffith,M., Griffith,O.L.,\n            Krzywinski,M.I., Liao,N., Morin,R., Palmquist,D., Petrescu,A.S.,\n            Skalska,U., Smailus,D.E., Stott,J.M., Schnerch,A., Schein,J.E.,\n            Jones,S.J., Holt,R.A., Baross,A., Marra,M.A., Clifton,S.,\n            Makowski,K.A., Bosak,S. and Malek,J.\n  CONSRTM   MGC Project Team\n  TITLE     The status, quality, and expansion of the NIH full-length cDNA\n            project: the Mammalian Gene Collection (MGC)\n  JOURNAL   Genome Res 14 (10B), 2121-2127 (2004)\n   PUBMED   15489334\n  REMARK    NUCLEOTIDE SEQUENCE [LARGE SCALE MRNA] (ISOFORM 1).;\n            STRAIN=Czech II; TISSUE=Eye, and Mammary tumor\n            Erratum:[Genome Res. 2006 Jun;16(6):804. Morrin, Ryan [corrected to\n            Morin, Ryan]]\nREFERENCE   6  (residues 1 to 410)\n  AUTHORS   Lubec,G., Klug,S., Sunyer,B. and Chen,W.-Q.\n  TITLE     Direct Submission\n  JOURNAL   Submitted (??-JAN-2009) to UniProtKB\n  REMARK    PROTEIN SEQUENCE OF 133-144 AND 375-405, AND IDENTIFICATION BY MASS\n            SPECTROMETRY.;\n            STRAIN=OF1; TISSUE=Hippocampus\nREFERENCE   7  (residues 1 to 410)\n  AUTHORS   Caspi,M., Atlas,R., Kantor,A., Sapir,T. and Reiner,O.\n  TITLE     Interaction between LIS1 and doublecortin, two lissencephaly gene\n            products\n  JOURNAL   Hum Mol Genet 9 (15), 2205-2213 (2000)\n   PUBMED   11001923\n  REMARK    INTERACTION WITH DCX, AND DEVELOPMENTAL STAGE.\nREFERENCE   8  (residues 1 to 410)\n  AUTHORS   Smith,D.S., Niethammer,M., Ayala,R., Zhou,Y., Gambello,M.J.,\n            Wynshaw-Boris,A. and Tsai,L.H.\n  TITLE     Regulation of cytoplasmic dynein behaviour and microtubule\n            organization by mammalian Lis1\n  JOURNAL   Nat Cell Biol 2 (11), 767-775 (2000)\n   PUBMED   11056530\n  REMARK    FUNCTION, INTERACTION WITH DYNEIN AND DYNACTIN, SUBCELLULAR\n            LOCATION, AND DEVELOPMENTAL STAGE.\nREFERENCE   9  (residues 1 to 410)\n  AUTHORS   Feng,Y., Olson,E.C., Stukenberg,P.T., Flanagan,L.A., Kirschner,M.W.\n            and Walsh,C.A.\n  TITLE     LIS1 regulates CNS lamination by interacting with mNudE, a central\n            component of the centrosome\n  JOURNAL   Neuron 28 (3), 665-679 (2000)\n   PUBMED   11163258\n  REMARK    INTERACTION WITH NDE1.\nREFERENCE   10 (residues 1 to 410)\n  AUTHORS   Sasaki,S., Shionoya,A., Ishida,M., Gambello,M.J., Yingling,J.,\n            Wynshaw-Boris,A. and Hirotsune,S.\n  TITLE     A LIS1/NUDEL/cytoplasmic dynein heavy chain complex in the\n            developing and adult nervous system\n  JOURNAL   Neuron 28 (3), 681-696 (2000)\n   PUBMED   11163259\n  REMARK    SELF-ASSOCIATION, INTERACTION WITH DYNEIN AND NDEL1, SUBCELLULAR\n            LOCATION, TISSUE SPECIFICITY, DEVELOPMENTAL STAGE, AND MUTAGENESIS\n            OF HIS-149.\nREFERENCE   11 (residues 1 to 410)\n  AUTHORS   Niethammer,M., Smith,D.S., Ayala,R., Peng,J., Ko,J., Lee,M.S.,\n            Morabito,M. and Tsai,L.H.\n  TITLE     NUDEL is a novel Cdk5 substrate that associates with LIS1 and\n            cytoplasmic dynein\n  JOURNAL   Neuron 28 (3), 697-711 (2000)\n   PUBMED   11163260\n  REMARK    INTERACTION WITH NDEL1, SUBCELLULAR LOCATION, AND TISSUE\n            SPECIFICITY.\nREFERENCE   12 (residues 1 to 410)\n  AUTHORS   Sweeney,K.J., Prokscha,A. and Eichele,G.\n  TITLE     NudE-L, a novel Lis1-interacting protein, belongs to a family of\n            vertebrate coiled-coil proteins\n  JOURNAL   Mech Dev 101 (1-2), 21-33 (2001)\n   PUBMED   11231056\n  REMARK    INTERACTION WITH NDEL1, DEVELOPMENTAL STAGE, AND MUTAGENESIS OF\n            HIS-149; SER-152 AND SER-169.;\n            STRAIN=Swiss Webster / NIH\nREFERENCE   13 (residues 1 to 410)\n  AUTHORS   Cahana,A., Escamez,T., Nowakowski,R.S., Hayes,N.L., Giacobini,M.,\n            von Holst,A., Shmueli,O., Sapir,T., McConnell,S.K., Wurst,W.,\n            Martinez,S. and Reiner,O.\n  TITLE     Targeted mutagenesis of Lis1 disrupts cortical development and LIS1\n            homodimerization\n  JOURNAL   Proc Natl Acad Sci U S A 98 (11), 6429-6434 (2001)\n   PUBMED   11344260\n  REMARK    FUNCTION, SELF-ASSOCIATION, INTERACTION WITH PAFAH1B2 AND PAFAH1B3,\n            AND DEVELOPMENTAL STAGE.\nREFERENCE   14 (residues 1 to 410)\n  AUTHORS   Coquelle,F.M., Caspi,M., Cordelieres,F.P., Dompierre,J.P.,\n            Dujardin,D.L., Koifman,C., Martin,P., Hoogenraad,C.C.,\n            Akhmanova,A., Galjart,N., De Mey,J.R. and Reiner,O.\n  TITLE     LIS1, CLIP-170's key to the dynein/dynactin pathway\n  JOURNAL   Mol Cell Biol 22 (9), 3089-3102 (2002)\n   PUBMED   11940666\n  REMARK    SUBCELLULAR LOCATION, AND INTERACTION WITH RSN.\nREFERENCE   15 (residues 1 to 410)\n  AUTHORS   Tokuoka,S.M., Ishii,S., Kawamura,N., Satoh,M., Shimada,A.,\n            Sasaki,S., Hirotsune,S., Wynshaw-Boris,A. and Shimizu,T.\n  TITLE     Involvement of platelet-activating factor and LIS1 in neuronal\n            migration\n  JOURNAL   Eur J Neurosci 18 (3), 563-570 (2003)\n   PUBMED   12911752\n  REMARK    FUNCTION.\nREFERENCE   16 (residues 1 to 410)\n  AUTHORS   Koizumi,H., Yamaguchi,N., Hattori,M., Ishikawa,T.O., Aoki,J.,\n            Taketo,M.M., Inoue,K. and Arai,H.\n  TITLE     Targeted disruption of intracellular type I platelet activating\n            factor-acetylhydrolase catalytic subunits causes severe impairment\n            in spermatogenesis\n  JOURNAL   J Biol Chem 278 (14), 12489-12494 (2003)\n   PUBMED   12551946\n  REMARK    SUBCELLULAR LOCATION, AND TISSUE SPECIFICITY.\nREFERENCE   17 (residues 1 to 410)\n  AUTHORS   Caspi,M., Coquelle,F.M., Koifman,C., Levy,T., Arai,H., Aoki,J., De\n            Mey,J.R. and Reiner,O.\n  TITLE     LIS1 missense mutations: variable phenotypes result from\n            unpredictable alterations in biochemical and cellular properties\n  JOURNAL   J Biol Chem 278 (40), 38740-38748 (2003)\n   PUBMED   12885786\n  REMARK    SELF-ASSOCIATION, INTERACTION WITH NUDC; NDE1; NDEL1; PAFAH1B2;\n            PAFAH1B3 AND RSN, SUBCELLULAR LOCATION, AND MUTAGENESIS OF PHE-31;\n            HIS-149; GLY-162; SER-169 AND ASP-317.\nREFERENCE   18 (residues 1 to 410)\n  AUTHORS   Dujardin,D.L., Barnhart,L.E., Stehman,S.A., Gomes,E.R.,\n            Gundersen,G.G. and Vallee,R.B.\n  TITLE     A role for cytoplasmic dynein and LIS1 in directed cell movement\n  JOURNAL   J Cell Biol 163 (6), 1205-1211 (2003)\n   PUBMED   14691133\n  REMARK    FUNCTION, AND SUBCELLULAR LOCATION.\nREFERENCE   19 (residues 1 to 410)\n  AUTHORS   Kholmanskikh,S.S., Dobrin,J.S., Wynshaw-Boris,A., Letourneau,P.C.\n            and Ross,M.E.\n  TITLE     Disregulated RhoGTPases and actin cytoskeleton contribute to the\n            migration defect in Lis1-deficient neurons\n  JOURNAL   J Neurosci 23 (25), 8673-8681 (2003)\n   PUBMED   14507966\n  REMARK    FUNCTION.\nREFERENCE   20 (residues 1 to 410)\n  AUTHORS   Cahana,A., Jin,X.L., Reiner,O., Wynshaw-Boris,A. and O'Neill,C.\n  TITLE     A study of the nature of embryonic lethality in LIS1-/- mice\n  JOURNAL   Mol Reprod Dev 66 (2), 134-142 (2003)\n   PUBMED   12950100\n  REMARK    SUBCELLULAR LOCATION, AND DEVELOPMENTAL STAGE.\nREFERENCE   21 (residues 1 to 410)\n  AUTHORS   Toyo-oka,K., Shionoya,A., Gambello,M.J., Cardoso,C., Leventer,R.,\n            Ward,H.L., Ayala,R., Tsai,L.H., Dobyns,W., Ledbetter,D.,\n            Hirotsune,S. and Wynshaw-Boris,A.\n  TITLE     14-3-3epsilon is important for neuronal migration by binding to\n            NUDEL: a molecular explanation for Miller-Dieker syndrome\n  JOURNAL   Nat Genet 34 (3), 274-285 (2003)\n   PUBMED   12796778\n  REMARK    FUNCTION, AND SUBCELLULAR LOCATION.\nREFERENCE   22 (residues 1 to 410)\n  AUTHORS   Assadi,A.H., Zhang,G., Beffert,U., McNeil,R.S., Renfro,A.L.,\n            Niu,S., Quattrocchi,C.C., Antalffy,B.A., Sheldon,M.,\n            Armstrong,D.D., Wynshaw-Boris,A., Herz,J., D'Arcangelo,G. and\n            Clark,G.D.\n  TITLE     Interaction of reelin signaling and Lis1 in brain development\n  JOURNAL   Nat Genet 35 (3), 270-276 (2003)\n   PUBMED   14578885\n  REMARK    FUNCTION, INTERACTION WITH DAB1, SUBCELLULAR LOCATION, AND\n            MUTAGENESIS OF HIS-149.\nREFERENCE   23 (residues 1 to 410)\n  AUTHORS   Yamaguchi,N., Takanezawa,Y., Koizumi,H., Umezu-Goto,M., Aoki,J. and\n            Arai,H.\n  TITLE     Expression of NUDEL in manchette and its implication in\n            spermatogenesis\n  JOURNAL   FEBS Lett 566 (1-3), 71-76 (2004)\n   PUBMED   15147871\n  REMARK    SUBCELLULAR LOCATION, AND DEVELOPMENTAL STAGE.\nREFERENCE   24 (residues 1 to 410)\n  AUTHORS   Tanaka,T., Serneo,F.F., Higgins,C., Gambello,M.J., Wynshaw-Boris,A.\n            and Gleeson,J.G.\n  TITLE     Lis1 and doublecortin function with dynein to mediate coupling of\n            the nucleus to the centrosome in neuronal migration\n  JOURNAL   J Cell Biol 165 (5), 709-721 (2004)\n   PUBMED   15173193\n  REMARK    FUNCTION, INTERACTION WITH DYNEIN AND DCX, AND SUBCELLULAR\n            LOCATION.\nREFERENCE   25 (residues 1 to 410)\n  AUTHORS   Shu,T., Ayala,R., Nguyen,M.D., Xie,Z., Gleeson,J.G. and Tsai,L.H.\n  TITLE     Ndel1 operates in a common pathway with LIS1 and cytoplasmic dynein\n            to regulate cortical neuronal positioning\n  JOURNAL   Neuron 44 (2), 263-277 (2004)\n   PUBMED   15473966\n  REMARK    FUNCTION, INTERACTION WITH DYNEIN AND DCX, SUBCELLULAR LOCATION,\n            AND DEVELOPMENTAL STAGE.\nREFERENCE   26 (residues 1 to 410)\n  AUTHORS   Feng,Y. and Walsh,C.A.\n  TITLE     Mitotic spindle regulation by Nde1 controls cerebral cortical size\n  JOURNAL   Neuron 44 (2), 279-293 (2004)\n   PUBMED   15473967\n  REMARK    INTERACTION WITH NDE1.\nREFERENCE   27 (residues 1 to 410)\n  AUTHORS   Gerlitz,G., Darhin,E., Giorgio,G., Franco,B. and Reiner,O.\n  TITLE     Novel functional features of the Lis-H domain: role in protein\n            dimerization, half-life and cellular localization\n  JOURNAL   Cell Cycle 4 (11), 1632-1640 (2005)\n   PUBMED   16258276\n  REMARK    SELF-ASSOCIATION, DOMAIN, AND MUTAGENESIS OF ILE-15 AND LEU-19.\nREFERENCE   28 (residues 1 to 410)\n  AUTHORS   Toyo-Oka,K., Sasaki,S., Yano,Y., Mori,D., Kobayashi,T.,\n            Toyoshima,Y.Y., Tokuoka,S.M., Ishii,S., Shimizu,T., Muramatsu,M.,\n            Hiraiwa,N., Yoshiki,A., Wynshaw-Boris,A. and Hirotsune,S.\n  TITLE     Recruitment of katanin p60 by phosphorylated NDEL1, an LIS1\n            interacting protein, is essential for mitotic cell division and\n            neuronal migration\n  JOURNAL   Hum Mol Genet 14 (21), 3113-3128 (2005)\n   PUBMED   16203747\n  REMARK    FUNCTION, AND INTERACTION WITH KATNB1.\nREFERENCE   29 (residues 1 to 410)\n  AUTHORS   Sasaki,S., Mori,D., Toyo-oka,K., Chen,A., Garrett-Beal,L.,\n            Muramatsu,M., Miyagawa,S., Hiraiwa,N., Yoshiki,A., Wynshaw-Boris,A.\n            and Hirotsune,S.\n  TITLE     Complete loss of Ndel1 results in neuronal migration defects and\n            early embryonic lethality\n  JOURNAL   Mol Cell Biol 25 (17), 7812-7827 (2005)\n   PUBMED   16107726\n  REMARK    FUNCTION.\nREFERENCE   30 (residues 1 to 410)\n  AUTHORS   Soukoulis,V., Reddy,S., Pooley,R.D., Feng,Y., Walsh,C.A. and\n            Bader,D.M.\n  TITLE     Cytoplasmic LEK1 is a regulator of microtubule function through its\n            interaction with the LIS1 pathway\n  JOURNAL   Proc Natl Acad Sci U S A 102 (24), 8549-8554 (2005)\n   PUBMED   15939891\n  REMARK    SUBCELLULAR LOCATION.\nREFERENCE   31 (residues 1 to 410)\n  AUTHORS   Mesngon,M.T., Tarricone,C., Hebbar,S., Guillotte,A.M.,\n            Schmitt,E.W., Lanier,L., Musacchio,A., King,S.J. and Smith,D.S.\n  TITLE     Regulation of cytoplasmic dynein ATPase by Lis1\n  JOURNAL   J Neurosci 26 (7), 2132-2139 (2006)\n   PUBMED   16481446\n  REMARK    FUNCTION, AND INTERACTION WITH DYNEIN AND PAFAH1B2.\nREFERENCE   32 (residues 1 to 410)\n  AUTHORS   Guo,J., Yang,Z., Song,W., Chen,Q., Wang,F., Zhang,Q. and Zhu,X.\n  TITLE     Nudel contributes to microtubule anchoring at the mother centriole\n            and is involved in both dynein-dependent and -independent\n            centrosomal protein assembly\n  JOURNAL   Mol Biol Cell 17 (2), 680-689 (2006)\n   PUBMED   16291865\n  REMARK    INTERACTION WITH NDEL1, AND SUBCELLULAR LOCATION.\nREFERENCE   33 (residues 1 to 410)\n  AUTHORS   Kholmanskikh,S.S., Koeller,H.B., Wynshaw-Boris,A., Gomez,T.,\n            Letourneau,P.C. and Ross,M.E.\n  TITLE     Calcium-dependent interaction of Lis1 with IQGAP1 and Cdc42\n            promotes neuronal motility\n  JOURNAL   Nat Neurosci 9 (1), 50-57 (2006)\n   PUBMED   16369480\n  REMARK    FUNCTION, AND INTERACTION WITH IQGAP1 AND RSN.\nREFERENCE   34 (residues 1 to 410)\n  AUTHORS   Zhang,G., Assadi,A.H., McNeil,R.S., Beffert,U., Wynshaw-Boris,A.,\n            Herz,J., Clark,G.D. and D'Arcangelo,G.\n  TITLE     The Pafah1b complex interacts with the reelin receptor VLDLR\n  JOURNAL   PLoS One 2 (2), e252 (2007)\n   PUBMED   17330141\n  REMARK    DISRUPTION PHENOTYPE, AND FUNCTION.\n            Publication Status: Online-Only\nREFERENCE   35 (residues 1 to 410)\n  AUTHORS   Huttlin,E.L., Jedrychowski,M.P., Elias,J.E., Goswami,T., Rad,R.,\n            Beausoleil,S.A., Villen,J., Haas,W., Sowa,M.E. and Gygi,S.P.\n  TITLE     A tissue-specific atlas of mouse protein phosphorylation and\n            expression\n  JOURNAL   Cell 143 (7), 1174-1189 (2010)\n   PUBMED   21183079\n  REMARK    IDENTIFICATION BY MASS SPECTROMETRY [LARGE SCALE ANALYSIS].;\n            TISSUE=Brain, Brown adipose tissue, Heart, Kidney, Liver, Lung,\n            Pancreas, Spleen, and Testis\nREFERENCE   36 (residues 1 to 410)\n  AUTHORS   Jodoin,J.N., Shboul,M., Sitaram,P., Zein-Sabatto,H., Reversade,B.,\n            Lee,E. and Lee,L.A.\n  TITLE     Human Asunder promotes dynein recruitment and centrosomal tethering\n            to the nucleus at mitotic entry\n  JOURNAL   Mol Biol Cell 23 (24), 4713-4724 (2012)\n   PUBMED   23097494\n  REMARK    INTERACTION WITH INTS13.\nREFERENCE   37 (residues 1 to 410)\n  AUTHORS   Tarricone,C., Perrina,F., Monzani,S., Massimiliano,L., Kim,M.H.,\n            Derewenda,Z.S., Knapp,S., Tsai,L.H. and Musacchio,A.\n  TITLE     Coupling PAF signaling to dynein regulation: structure of LIS1 in\n            complex with PAF-acetylhydrolase\n  JOURNAL   Neuron 44 (5), 809-821 (2004)\n   PUBMED   15572112\n  REMARK    X-RAY CRYSTALLOGRAPHY (3.4 ANGSTROMS) IN COMPLEX WITH PAFAH1B2,\n            INTERACTION WITH NDEL1, AND MUTAGENESIS OF ARG-238.\nREFERENCE   38 (residues 1 to 410)\n  AUTHORS   Kim,M.H., Cooper,D.R., Oleksy,A., Devedjiev,Y., Derewenda,U.,\n            Reiner,O., Otlewski,J. and Derewenda,Z.S.\n  TITLE     The structure of the N-terminal domain of the product of the\n            lissencephaly gene Lis1 and its functional implications\n  JOURNAL   Structure 12 (6), 987-998 (2004)\n   PUBMED   15274919\n  REMARK    X-RAY CRYSTALLOGRAPHY (1.75 ANGSTROMS) OF 1-86, DOMAIN, AND\n            CIRCULAR DICHROISM ANALYSIS.\nCOMMENT     [FUNCTION] Regulatory subunit (beta subunit) of the cytosolic type\n            I platelet-activating factor (PAF) acetylhydrolase (PAF-AH (I)), an\n            enzyme that catalyzes the hydrolyze of the acetyl group at the sn-2\n            position of PAF and its analogs and participates in PAF\n            inactivation. Regulates the PAF-AH (I) activity in a catalytic\n            dimer composition-dependent manner (By similarity). Positively\n            regulates the activity of the minus-end directed microtubule motor\n            protein dynein. May enhance dynein-mediated microtubule sliding by\n            targeting dynein to the microtubule plus end. Required for several\n            dynein- and microtubule-dependent processes such as the maintenance\n            of Golgi integrity, the peripheral transport of microtubule\n            fragments and the coupling of the nucleus and centrosome. Required\n            during brain development for the proliferation of neuronal\n            precursors and the migration of newly formed neurons from the\n            ventricular/subventricular zone toward the cortical plate. Neuronal\n            migration involves a process called nucleokinesis, whereby\n            migrating cells extend an anterior process into which the nucleus\n            subsequently translocates. During nucleokinesis dynein at the\n            nuclear surface may translocate the nucleus towards the centrosome\n            by exerting force on centrosomal microtubules. Also required for\n            proper activation of Rho GTPases and actin polymerization at the\n            leading edge of locomoting cerebellar neurons and postmigratory\n            hippocampal neurons in response to calcium influx triggered via\n            NMDA receptors. May also play a role in other forms of cell\n            locomotion including the migration of fibroblasts during wound\n            healing. Non-catalytic subunit of an acetylhydrolase complex which\n            inactivates platelet-activating factor (PAF) by removing the acetyl\n            group at the SN-2 position. Required for dynein recruitment to\n            microtubule plus ends and BICD2-bound cargos (By similarity). May\n            modulate the Reelin pathway through interaction of the PAF-AH (I)\n            catalytic dimer with VLDLR (PubMed:17330141).\n            {ECO:0000250|UniProtKB:P43033, ECO:0000250|UniProtKB:P43034,\n            ECO:0000255|HAMAP-Rule:MF_03141, ECO:0000269|PubMed:11056530,\n            ECO:0000269|PubMed:11344260, ECO:0000269|PubMed:12796778,\n            ECO:0000269|PubMed:12911752, ECO:0000269|PubMed:14507966,\n            ECO:0000269|PubMed:14578885, ECO:0000269|PubMed:14691133,\n            ECO:0000269|PubMed:15173193, ECO:0000269|PubMed:15473966,\n            ECO:0000269|PubMed:16107726, ECO:0000269|PubMed:16203747,\n            ECO:0000269|PubMed:16369480, ECO:0000269|PubMed:16481446,\n            ECO:0000269|PubMed:17330141}.\n            [SUBUNIT] Component of the cytosolic PAF-AH (I) heterotetrameric\n            enzyme, which is composed of PAFAH1B1 (beta), PAFAH1B2 (alpha2) and\n            PAFAH1B3 (alpha1) subunits. The catalytic activity of the enzyme\n            resides in the alpha1 (PAFAH1B3) and alpha2 (PAFAH1B2) subunits,\n            whereas the beta subunit (PAFAH1B1) has regulatory activity. Trimer\n            formation is not essential for the catalytic activity. Interacts\n            with the catalytic dimer of PAF-AH (I) heterotetrameric enzyme:\n            interacts with PAFAH1B2 homodimer (alpha2/alpha2 homodimer),\n            PAFAH1B3 homodimer (alpha1/alpha1 homodimer) and PAFAH1B2-PAFAH1B3\n            heterodimer (alpha2/alpha1 heterodimer) (PubMed:11344260,\n            PubMed:12885786, PubMed:16481446). Can self-associate\n            (PubMed:11163259, PubMed:11344260, PubMed:12885786,\n            PubMed:15572112). Interacts with DCX, dynein, dynactin, IQGAP1,\n            KATNB1, NDE1, NDEL1, NUDC, and RSN (PubMed:11001923,\n            PubMed:11056530, PubMed:11163258, PubMed:11163259, PubMed:11163260,\n            PubMed:11231056, PubMed:11940666, PubMed:12885786, PubMed:15173193,\n            PubMed:15473966, PubMed:15473967, PubMed:15572112, PubMed:16203747,\n            PubMed:16291865, PubMed:16481446, PubMed:16369480). Interacts with\n            DAB1 when DAB1 is phosphorylated in response to RELN/reelin\n            signaling (PubMed:14578885). Interacts with INTS13\n            (PubMed:23097494). Interacts with DCDC1 (By similarity). Interacts\n            with DISC1, and this interaction is enhanced by NDEL1 (By\n            similarity). {ECO:0000250|UniProtKB:P43034,\n            ECO:0000269|PubMed:11001923, ECO:0000269|PubMed:11056530,\n            ECO:0000269|PubMed:11163258, ECO:0000269|PubMed:11163259,\n            ECO:0000269|PubMed:11163260, ECO:0000269|PubMed:11231056,\n            ECO:0000269|PubMed:11344260, ECO:0000269|PubMed:11940666,\n            ECO:0000269|PubMed:12885786, ECO:0000269|PubMed:14578885,\n            ECO:0000269|PubMed:15173193, ECO:0000269|PubMed:15473966,\n            ECO:0000269|PubMed:15473967, ECO:0000269|PubMed:15572112,\n            ECO:0000269|PubMed:16203747, ECO:0000269|PubMed:16291865,\n            ECO:0000269|PubMed:16369480, ECO:0000269|PubMed:16481446,\n            ECO:0000269|PubMed:23097494}.\n            [INTERACTION] P63005; Q9ERR1: Ndel1; NbExp=9; IntAct=EBI-917499,\n            EBI-646668; P63005; Q61205: Pafah1b3; NbExp=2; IntAct=EBI-917499,\n            EBI-1007637; P63005; Q9GZM8: NDEL1; Xeno; NbExp=4;\n            IntAct=EBI-917499, EBI-928842.\n            [SUBCELLULAR LOCATION] Cytoplasm, cytoskeleton. Cytoplasm,\n            cytoskeleton, microtubule organizing center, centrosome. Cytoplasm,\n            cytoskeleton, spindle. Nucleus membrane\n            {ECO:0000255|HAMAP-Rule:MF_03141}. Note=May localize to the nuclear\n            membrane (By similarity). Localizes to the plus end of microtubules\n            and to the centrosome. Redistributes to axons during neuronal\n            development. Also localizes to the microtubules of the manchette in\n            elongating spermatids and to the meiotic spindle in spermatocytes.\n            {ECO:0000250}.\n            [ALTERNATIVE PRODUCTS] Event=Alternative splicing; Named\n            isoforms=2; Name=1; IsoId=P63005-1, P43035-1; Sequence=Displayed;\n            Name=2; IsoId=P63005-2, P43035-2; Sequence=VSP_006778.\n            [TISSUE SPECIFICITY] Highly expressed in brain, particularly the\n            hippocampus and the olfactory bulb. Also highly expressed in\n            testis, including all seminiferous tubule cell types, all types of\n            spermatogenic and Sertoli cells, and meiotically dividing and\n            elongating spermatids. Expressed at lower levels in heart, kidney,\n            large intestine, liver, lung, ovary, small intestine and spleen.\n            {ECO:0000269|PubMed:11163259, ECO:0000269|PubMed:11163260,\n            ECO:0000269|PubMed:12551946}.\n            [DEVELOPMENTAL STAGE] Embryonic expression begins prior to the\n            blastocyst stage, when maternally expressed protein is depleted. By\n            10.5 dpc, expression is abundant in the developing central and\n            peripheral nervous systems. Major sites of expression include the\n            neuroepithelium of the fore-, mid-, and hindbrain, the spinal cord,\n            the dorsal root and the cranial ganglia. By 13.5 dpc, highly\n            expressed in neuroblasts as well as postmitotic neurons of the\n            cortical plate. After completion of neuronal migration expression\n            remains high in the cortex. Also expressed in the testis from P8.\n            {ECO:0000269|PubMed:11001923, ECO:0000269|PubMed:11056530,\n            ECO:0000269|PubMed:11163259, ECO:0000269|PubMed:11231056,\n            ECO:0000269|PubMed:11344260, ECO:0000269|PubMed:12950100,\n            ECO:0000269|PubMed:15147871, ECO:0000269|PubMed:15473966}.\n            [DOMAIN] Dimerization mediated by the LisH domain may be required\n            to activate dynein. {ECO:0000255|HAMAP-Rule:MF_03141,\n            ECO:0000269|PubMed:15274919, ECO:0000269|PubMed:16258276}.\n            [DISRUPTION PHENOTYPE] Double heterozygous PAFAH1B1 and homozygous\n            VLDLR knockout mice present no obvious cortical layering defects\n            (PubMed:17330141). Double heterozygous PAFAH1B1 and homozygous LRP8\n            knockout mice display a reeler-like phenotype in the forebrain\n            consisting of the inversion of cortical layers and hippocampal\n            disorganization (PubMed:17330141). {ECO:0000269|PubMed:17330141}.\n            [MISCELLANEOUS] Originally the subunits of the type I\n            platelet-activating factor (PAF) acetylhydrolase was named alpha\n            (PAFAH1B1), beta (PAFAH1B2) and gamma (PAFAH1B3) (By similarity).\n            Now these subunits have been renamed beta (PAFAH1B1), alpha2\n            (PAFAH1B2) and alpha1 (PAFAH1B3) respectively (By similarity).\n            {ECO:0000250|UniProtKB:P43034, ECO:0000250|UniProtKB:P68402,\n            ECO:0000250|UniProtKB:Q15102, ECO:0000250|UniProtKB:Q29460}.\n            [SIMILARITY] Belongs to the WD repeat LIS1/nudF family.\n            {ECO:0000255|HAMAP-Rule:MF_03141}.\nFEATURES             Location/Qualifiers\n     source          1..410\n                     /organism=\"Mus musculus\"\n                     /db_xref=\"taxon:10090\"\n     gene            1..410\n                     /gene=\"Pafah1b1\"\n                     /gene_synonym=\"Lis-1\"\n                     /gene_synonym=\"Lis1\"\n                     /gene_synonym=\"Pafaha\"\n     Protein         1..410\n                     /product=\"Platelet-activating factor acetylhydrolase IB\n                     subunit beta\"\n                     /note=\"Lissencephaly-1 protein; PAF acetylhydrolase 45 kDa\n                     subunit; PAF-AH alpha; LIS-1; PAF-AH 45 kDa subunit; PAFAH\n                     alpha\"\n                     /UniProtKB_evidence=\"Evidence at protein level\"\n     Region          1..410\n                     /region_name=\"Mature chain\"\n                     /note=\"Platelet-activating factor acetylhydrolase IB\n                     subunit beta. /id=PRO_0000051063.\"\n     Region          1..102\n                     /region_name=\"Region of interest in the sequence\"\n                     /note=\"Interaction with NDEL1.\"\n     Region          1..66\n                     /region_name=\"Region of interest in the sequence\"\n                     /note=\"Interaction with NDE1.\n                     /evidence=ECO:0000255|HAMAP-Rule:MF_03141.\"\n     Region          1..38\n                     /region_name=\"Region of interest in the sequence\"\n                     /note=\"Required for self-association and interaction with\n                     PAFAH1B2 and PAFAH1B3.\"\n     Region          5..21\n                     /region_name=\"Helical region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1UUJ.\"\n     Region          7..39\n                     /region_name=\"Domain\"\n                     /note=\"LisH. /evidence=ECO:0000255|HAMAP-Rule:MF_03141.\"\n     Region          10..34\n                     /region_name=\"LisH\"\n                     /note=\"pfam08513\"\n                     /db_xref=\"CDD:369916\"\n     Site            15\n                     /site_type=\"mutagenized\"\n                     /note=\"I->R: Impairs self-association and reduces protein\n                     stability. /evidence=ECO:0000269|PubMed:16258276.\"\n     Site            19\n                     /site_type=\"mutagenized\"\n                     /note=\"L->R: Impairs self-association, reduces protein\n                     stability and promotes localization to actin fibers.\n                     /evidence=ECO:0000269|PubMed:16258276.\"\n     Region          25..34\n                     /region_name=\"Helical region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1UUJ.\"\n     Site            31\n                     /site_type=\"mutagenized\"\n                     /note=\"F->S: Abrogates interaction with NDLE1, NUDC,\n                     PAFAH1B3 and RSN. Also impairs localization to centrosomes\n                     and microtubule plus ends. Reduces protein stability.\n                     /evidence=ECO:0000269|PubMed:12885786.\"\n     Region          41..47\n                     /region_name=\"Helical region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1UUJ.\"\n     Region          50..55\n                     /region_name=\"Helical region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1UUJ.\"\n     Site            53\n                     /site_type=\"acetylation\"\n                     /note=\"N6-acetyllysine.\n                     /evidence=ECO:0000250|UniProtKB:P43034.\"\n     Region          56..82\n                     /region_name=\"Coiled-coil region\"\n                     /note=\"/evidence=ECO:0000255|HAMAP-Rule:MF_03141.\"\n     Region          58..74\n                     /region_name=\"Helical region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1UUJ.\"\n     Region          83..410\n                     /region_name=\"Region of interest in the sequence\"\n                     /note=\"Interaction with dynein and dynactin.\n                     /evidence=ECO:0000255|HAMAP-Rule:MF_03141.\"\n     Region          103..105\n                     /region_name=\"Beta-strand region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Region          104..408\n                     /region_name=\"WD40\"\n                     /note=\"WD40 domain, found in a number of eukaryotic\n                     proteins that cover a wide variety of functions including\n                     adaptor/regulatory modules in signal transduction,\n                     pre-mRNA processing and cytoskeleton assembly; typically\n                     contains a GH dipeptide 11-24 residues from...; cd00200\"\n                     /db_xref=\"CDD:238121\"\n     Region          106..147\n                     /region_name=\"Repetitive region\"\n                     /note=\"WD 1.\"\n     Site            order(107,125,129,135..136,148..149,167,171,177..178,\n                     190..191,208,213,219..220,233,250,255,261..262,274..275,\n                     313,317,323..324,336..337,354,359,365..366,378..379,397,\n                     401,407..408)\n                     /site_type=\"active\"\n                     /note=\"structural tetrad [active]\"\n                     /db_xref=\"CDD:238121\"\n     Site            109\n                     /site_type=\"phosphorylation\"\n                     /note=\"Phosphoserine.\n                     /evidence=ECO:0000250|UniProtKB:P43034.\"\n     Region          111..148\n                     /region_name=\"WD40 repeat\"\n                     /note=\"WD40 repeat [structural motif]\"\n                     /db_xref=\"CDD:293791\"\n     Region          111..116\n                     /region_name=\"Beta-strand region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Region          118..130\n                     /region_name=\"Beta-strand region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Region          132..136\n                     /region_name=\"Beta-strand region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Region          137..139\n                     /region_name=\"Hydrogen bonded turn\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Region          144..146\n                     /region_name=\"Beta-strand region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Region          148..187\n                     /region_name=\"Repetitive region\"\n                     /note=\"WD 2.\"\n     Site            149\n                     /site_type=\"mutagenized\"\n                     /note=\"H->R: Abrogates self-association and interaction\n                     with DAB1, dynein, NDEL1, PAFAH1B3 and RSN. Also impairs\n                     localization to centrosomes and microtubule plus ends.\n                     Reduces protein stability.\n                     /evidence=ECO:0000269|PubMed:11163259,\n                     ECO:0000269|PubMed:11231056, ECO:0000269|PubMed:12885786,\n                     ECO:0000269|PubMed:14578885.\"\n     Site            152\n                     /site_type=\"mutagenized\"\n                     /note=\"S->W: Abrogates interaction with NDEL1.\n                     /evidence=ECO:0000269|PubMed:11231056.\"\n     Region          153..158\n                     /region_name=\"Beta-strand region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Region          154..190\n                     /region_name=\"WD40 repeat\"\n                     /note=\"WD40 repeat [structural motif]\"\n                     /db_xref=\"CDD:293791\"\n     Region          162..169\n                     /region_name=\"Beta-strand region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Site            162\n                     /site_type=\"mutagenized\"\n                     /note=\"G->S: Abrogates interaction with PAFAH1B3 and RSN.\n                     Also impairs localization to centrosomes and microtubule\n                     plus ends. Reduces protein stability.\n                     /evidence=ECO:0000269|PubMed:12885786.\"\n     Site            169\n                     /site_type=\"mutagenized\"\n                     /note=\"S->P: Abrogates interaction with NDEL1, PAFAH1B3\n                     and RSN. Also impairs localization to centrosomes and\n                     microtubule plus ends. Reduces protein stability.\n                     /evidence=ECO:0000269|PubMed:11231056,\n                     ECO:0000269|PubMed:12885786.\"\n     Region          176..178\n                     /region_name=\"Beta-strand region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Region          184..186\n                     /region_name=\"Beta-strand region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Region          190..229\n                     /region_name=\"Repetitive region\"\n                     /note=\"WD 3.\"\n     Region          195..231\n                     /region_name=\"WD40 repeat\"\n                     /note=\"WD40 repeat [structural motif]\"\n                     /db_xref=\"CDD:293791\"\n     Region          195..200\n                     /region_name=\"Beta-strand region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Region          202..211\n                     /region_name=\"Beta-strand region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Region          214..220\n                     /region_name=\"Beta-strand region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Region          221..223\n                     /region_name=\"Hydrogen bonded turn\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Region          226..231\n                     /region_name=\"Beta-strand region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Region          232..271\n                     /region_name=\"Repetitive region\"\n                     /note=\"WD 4.\"\n     Region          237..242\n                     /region_name=\"Beta-strand region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Region          238..273\n                     /region_name=\"WD40 repeat\"\n                     /note=\"WD40 repeat [structural motif]\"\n                     /db_xref=\"CDD:293791\"\n     Site            238\n                     /site_type=\"mutagenized\"\n                     /note=\"R->A: Abrogates interaction with PAFAH1B2.\n                     /evidence=ECO:0000269|PubMed:15572112.\"\n     Region          246..253\n                     /region_name=\"Beta-strand region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Region          258..262\n                     /region_name=\"Beta-strand region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Region          263..265\n                     /region_name=\"Hydrogen bonded turn\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Region          268..272\n                     /region_name=\"Beta-strand region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Region          274..333\n                     /region_name=\"Repetitive region\"\n                     /note=\"WD 5.\"\n     Region          279..335\n                     /region_name=\"WD40 repeat\"\n                     /note=\"WD40 repeat [structural motif]\"\n                     /db_xref=\"CDD:293791\"\n     Region          279..284\n                     /region_name=\"Beta-strand region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Region          289..295\n                     /region_name=\"Helical region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Region          310..315\n                     /region_name=\"Beta-strand region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Site            317\n                     /site_type=\"mutagenized\"\n                     /note=\"D->H: Abrogates interaction with NDEL1, PAFAH1B3\n                     and RSN. Also impairs localization to centrosomes and\n                     microtubule plus ends. Reduces protein stability.\n                     /evidence=ECO:0000269|PubMed:12885786.\"\n     Region          318..324\n                     /region_name=\"Beta-strand region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Region          325..328\n                     /region_name=\"Hydrogen bonded turn\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Region          329..335\n                     /region_name=\"Beta-strand region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Region          336..377\n                     /region_name=\"Repetitive region\"\n                     /note=\"WD 6.\"\n     Region          341..377\n                     /region_name=\"WD40 repeat\"\n                     /note=\"WD40 repeat [structural motif]\"\n                     /db_xref=\"CDD:293791\"\n     Region          341..346\n                     /region_name=\"Beta-strand region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Region          348..351\n                     /region_name=\"Beta-strand region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Region          353..357\n                     /region_name=\"Beta-strand region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Region          358..360\n                     /region_name=\"Hydrogen bonded turn\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Region          361..365\n                     /region_name=\"Beta-strand region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Region          367..409\n                     /region_name=\"Region of interest in the sequence\"\n                     /note=\"Interaction with DCX.\n                     /evidence=ECO:0000255|HAMAP-Rule:MF_03141.\"\n     Region          374..377\n                     /region_name=\"Beta-strand region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Region          378..410\n                     /region_name=\"Repetitive region\"\n                     /note=\"WD 7.\"\n     Region          383..407\n                     /region_name=\"WD40 repeat\"\n                     /note=\"WD40 repeat [structural motif]\"\n                     /db_xref=\"CDD:293791\"\n     Region          383..388\n                     /region_name=\"Beta-strand region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Region          387..410\n                     /region_name=\"Splicing variant\"\n                     /note=\"DFHKTAPYVVTGSVDQTVKVWECR -> GMYTL (in isoform 2).\n                     /evidence=ECO:0000303|PubMed:7751941. /id=VSP_006778.\"\n     Region          388..410\n                     /region_name=\"Region of interest in the sequence\"\n                     /note=\"Interaction with NDEL1.\"\n     Region          390..393\n                     /region_name=\"Beta-strand region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Region          395..399\n                     /region_name=\"Beta-strand region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\n     Region          402..407\n                     /region_name=\"Beta-strand region\"\n                     /note=\"/evidence=ECO:0007829|PDB:1VYH.\"\nORIGIN      \n        1 mvlsqrqrde lnraiadylr sngyeeaysv fkkeaeldmn eeldkkyagl lekkwtsvir\n       61 lqkkvmeles klneakeeft sggplgqkrd pkewiprppe kyalsghrsp vtrvifhpvf\n      121 svmvsaseda tikvwdyetg dfertlkght dsvqdisfdh sgkllascsa dmtiklwdfq\n      181 gfecirtmhg hdhnvssvai mpngdhivsa srdktikmwe vqtgycvktf tghrewvrmv\n      241 rpnqdgtlia scsndqtvrv wvvatkecka elrehehvve ciswapessy ssiseatgse\n      301 tkksgkpgpf llsgsrdkti kmwdvstgmc lmtlvghdnw vrgvlfhsgg kfilscaddk\n      361 tlrvwdyknk rcmktlnahe hfvtsldfhk tapyvvtgsv dqtvkvwecr\n//\n\n"

  ncbi.out <- list();
  for(p in missing.uniprots){
      printf("--- %s", p)
      p.trimmed <- sub("-.*$", "", p)  # "Q8BHN1-3" is isoform 3 of this protein.
      ncbi.out[[p]] <- entrez_fetch(db="protein", id=p.trimmed, rettype="text")
      }

  tbls <- list()
  for(i in seq_len(length(ncbi.out))){
     uniprot.id <- names(ncbi.out)[i]
     text.raw <- substr(ncbi.out[[i]], 1, 30)
     text.cooked.0 <- sub("LOCUS *", "", text.raw)
     text.cooked.1 <- sub(" *$", "", text.cooked.0)
     tokens <- strsplit(text.cooked.1, "_")[[1]]
     geneSymbol <- tokens[1]
     species <- tokens[2]
     tbls[[i]] <- data.frame(uniprot=uniprot.id, geneSymbol=geneSymbol, species=species)
     } # for i

  tbl.ncbi <- do.call(rbind, tbls)
  save(tbl.ncbi, file="missing.uniprots.looked.up.at.entrez.RData")

} # lookup.unmapped.uniprotids
#------------------------------------------------------------------------------------------------------------------------
merge.biomart.and.entrez.id.maps <- function()
{
   print(load("tbls.after.cleanup.RData"))
   print(load("tbl.uniprots.mapped.to.geneSymbol.RData"))
   print(load("missing.uniprots.looked.up.at.entrez.RData"))

   tbl.mart.ready <- tbl.mart[, c("uniprotswissprot", "hgnc_symbol")]
   colnames(tbl.mart.ready) <- c("uniprot", "geneSymbol")
   tbl.mart.ready$species <- "HUMAN"

   intersect(tbl.ncbi$uniprot, tbl.mart.ready$uniprot)  # nothing
   tbl.ids <- unique(rbind(tbl.ncbi, tbl.mart.ready))
   dim(tbl.ids)
     # we have them all:
   checkTrue(length(setdiff(ids, tbl.ids$uniprot), 0))

     # the id map need only include the proteins observed in ne1, ne2, and cyto
   tbl.ids <- subset(tbl.ids, uniprot %in% ids)
   dim(tbl.ids)

   save(tbl.ids, file="id.mapping.table.RData")

} # merge.biomart.and.entrez.id.maps
#------------------------------------------------------------------------------------------------------------------------
add.geneSymbols.to.protein.matrices <- function()
{
   print(load("id.mapping.table.RData"))
   print(load("tbls.after.cleanup.RData")) # "tbl.cyto" "tbl.ne1"  "tbl.ne2"

   checkTrue(all(tbl.cyto$Protein %in% tbl.ids$uniprot))
   checkTrue(all(tbl.ne1$Protein %in% tbl.ids$uniprot))
   checkTrue(all(tbl.ne2$Protein %in% tbl.ids$uniprot))

       # tbl.cyto first

   indices <- match(tbl.cyto$Protein, tbl.ids$uniprot)
   tbl.cyto$gene <- tbl.ids$geneSymbol[indices]
   tbl.cyto$species <- tbl.ids$species[indices]
   table(tbl.cyto$species)  # BOVIN CANLF HUMAN MOUSE   PIG   RAT
                            #     1     1  4649    45     1    25
   tbl.cyto <- subset(tbl.cyto, species=="HUMAN")

       # now tbl.ne1
   indices <- match(tbl.ne1$Protein, tbl.ids$uniprot)
   tbl.ne1$gene <- tbl.ids$geneSymbol[indices]
   tbl.ne1$species <- tbl.ids$species[indices]
   table(tbl.ne1$species) # CANLF HUMAN MOUSE   PIG   RAT
                          #     1  3898    38     2    19
   tbl.ne1 <- subset(tbl.ne1, species=="HUMAN")

       # now tbl.ne2
   indices <- match(tbl.ne2$Protein, tbl.ids$uniprot)
   tbl.ne2$gene <- tbl.ids$geneSymbol[indices]
   tbl.ne2$species <- tbl.ids$species[indices]
   table(tbl.ne2$species) # CANLF HUMAN MOUSE   PIG   RAT
                          #     1  2945    31     1    12
   tbl.ne2 <- subset(tbl.ne2, species=="HUMAN")

   dim(tbl.cyto)  # 4649   13
   dim(tbl.ne1)   # 3898   13
   dim(tbl.ne2)   # 2945   13
   save(tbl.cyto, tbl.ne1, tbl.ne2,
        file="cyto-ne1-ne2-trimmed-colnames-human-only-with-geneSymbols.RData")

} # add.geneSymbols.to.protein.matrices
#------------------------------------------------------------------------------------------------------------------------
venn.diagram <- function()
{
   require(VennDiagram)
   print(load("cyto-ne1-ne2-trimmed-colnames-human-only-with-geneSymbols.RData"))

   venn.diagram(list(cyto=tbl.cyto$gene, ne1=tbl.ne1$gene, ne2=tbl.ne2$gene),
                fill=c("red", "green", "blue"), filename="venn.tiff")

   cyto <- tbl.cyto$gene
   ne1  <- tbl.ne1$gene
   ne2  <- tbl.ne2$gene

   all.unique <- length(unique(c(cyto, ne1, ne2)))
   cyto.unique <- length(setdiff(cyto, c(ne1, ne2)))/all.unique
   ne1.unique  <- length(setdiff(ne1, c(cyto, ne2)))/all.unique
   ne2.unique  <- length(setdiff(ne2, c(cyto, ne1)))/all.unique

   cyto.ne1.ne2.shared <- length(intersect(cyto, intersect(ne1, ne2)))/all.unique
   cyto.ne1.shared <- length(intersect(cyto, ne1))/all.unique
   cyto.ne2.shared <- length(intersect(cyto, ne2))/all.unique

   ne1.ne2.shared <- length(intersect(ne1, ne2))/all.unique

   v <- venneuler(c(A=cyto.unique, B=ne1.unique, C=ne2.unique, "A&B"=cyto.ne1.shared, "A&C"=cyto.ne2.shared, "B&C"=ne1.ne2.shared, "A&B&C"=cyto.ne1.ne2.shared))
   plot(v)

   m <- data.frame(cyto=tbl.cyto$gene,  ne1=tbl.ne1$gene, ne2=tbl.ne2$gene)
   v <- venneuler(m > 0)
   plot(v)

} # venn.diagram
#------------------------------------------------------------------------------------------------------------------------

# missing.map <- list(   # from manual search of https://www.ncbi.nlm.nih.gov/protein/?
#     "P11884" = "ALDH", # rat
#     "Q58FG0" = "HSP90AA5P", # human
#     "P06685" = "AT1A1", # rat
#     "O88600" = "HSP74", # rat
#     "Q86V81" = "THOC4", # human
#     "P47758" = "SRPRB", # mouse
#     "P49256" = "LMAN2", # wolf
#     "Q58FF6" = "H90B4", # human putative hsp 90 beta 4
#     "P63005" = "LIS1", # MOUSE
#     "P35433" = "", #
#     "Q14568" = "", #
#     "Q9QUL6" = "", #
#     "Q61576" = "", #
#     "Q61879" = "", #
#     "Q63151" = "", #
#     "P01834" = "", #
#     "P01137" = "", #
#     "Q08369" = "", #
#     "Q02248" = "", #
#     "Q62018" = "", #
#     "Q9WU70" = "", #
#     "Q9Z0P5" = "", #
#     "P24928" = "", #
#     "P00388" = "", #
#     "P19156" = "", #
#     "Q58FF8" = "", #
#     "P47757" = "", #
#     "Q9JMG7-2" = "", #
#     "O35547" = "", #
#     "Q08877" = "", #
#     "Q8C0E2" = "", #
#     "Q9R0P4" = "", #
#     "Q69ZR2" = "", #
#     "Q9D1P4" = "", #
#     "Q9JHJ0" = "", #
#     "Q9BRK5" = "", #
#     "Q62991" = "", #
#     "A2AF47" = "", #
#     "D2XV59" = "", #
#     "Q9BZK3" = "", #
#     "Q62824" = "", #
#     "P17427" = "", #
#     "C4AMC7" = "", #
#     "Q8BHS6" = "", #
#     "Q8BH64" = "", #
#     "Q8R3H7" = "", #
#     "Q3UH60" = "", #
#     "P35290" = "", #
#     "Q80W54" = "", #
#     "Q61171" = "", #
#     "Q8BX90" = "", #
#     "Q86UK7" = "", #
#     "P21575" = "", #
#     "P22091" = "", #
#     "Q9H270" = "", #
#     "Q8BZ47" = "", #
#     "Q02038" = "", #
#     "Q6A028" = "", #
#     "Q8BID6" = "", #
#     "Q8VIJ5" = "", #
#     "Q8BHN1-3" = "", #
#     "Q8VEJ4" = "", #
#     "Q62768" = "", #
#     "Q9QZS3" = "", #
#     "Q9Z247" = "", #
#     "Q64127" = "", #
#     "P28660" = "", #
#     "Q5T3I0" = "", #
#     "O35274" = "", #
#     "P70587" = "", #
#     "Q9CZV8" = "", #
#     "P46660" = "", #
#     "Q62771-2" = "", #
#     "Q80V26" = "", #
#     "Q96GQ7" = "", #
#     "Q9Z103" = "", #
#     "O70546" = "", #
#     "Q8NI36" = "", #
#     "Q6DN03" = "", #
#     "Q99104" = "", #
#     "Q8N954" = "", #
#     "Q63357" = "", #
#     "P58501" = "", #
#     "Q9CR47" = "", #
#     "Q80T69" = "", #
#     "Q9EQJ4" = "", #
#     "Q9UKL3" = "", #
#     "Q6PCM1" = "", #
#     "Q7TQN3" = "", #
#     "P17012" = "", #
#     "P51576" = "", #
#     "P67999" = "", #
#     "Q6IE36" = "", #
#     "A4FUI1" = "", #
#     "P97798-2" = "", #
#     "P35520-2" = "", #
#     "P28037" = "", #
#     "Q60954" = "", #
#     "P49185" = "", #
#     "P49615" = "", #
#     "Q8K2C8" = "", #
#     "Q8C4V4" = "", #
#     "O35095" = "", #
#     "Q7TSS2" = "", #
#     "P97924" = "", #
#     "P47728" = "", #
#     "Q8R1R3" = "", #
#     "Q96EY9" = "", #
#     "Q9HB19" = "", #
#     "Q8VEE1" = "", #
#     "Q6PAV2" = "", #
#     "Q5SSL4" = "", #
#     "Q6PFY1-2" = "", #
#     "Q9JME2" = "", #
#     "P70451" = "", #
#     "A8MPP1" = "", #
#     "Q96JP2" = "", #
#     "Q6PDX6" = "", #
#     "Q02563" = "", #
#     "Q8C1Y8" = "" #
#     )
#
# library(EnsDb.Hsapiens.v86)
# tbl.ensdb <- select(EnsDb.Hsapiens.v86, keys=ids, columns=c("UNIPROTID", "GENEID","SYMBOL"), keytype="UNIPROTID")
#
# all(ids %in%
#
# library(AnnotationHub)
# ah <- AnnotationHub()
# query(ah, pattern = c("Homo Sapiens", "EnsDb"))
#   #   AH89426  | Ensembl 103 EnsDb for Homo sapiens
#   #   AH95744  | Ensembl 104 EnsDb for Homo sapiens
#   #   AH98047  | Ensembl 105 EnsDb for Homo sapiens
#   #   AH100643 | Ensembl 106 EnsDb for Homo sapiens
#   #   AH104864 | Ensembl 107 EnsDb for Homo sapiens
#
# db <- ah[["AH104864"]]
# length(ids)      # 5790
# head(ids)
#
# > columns(db)
# #  [1] "CANONICALTRANSCRIPT" "DESCRIPTION"         "ENTREZID"            "EXONID"
# #  [5] "EXONIDX"             "EXONSEQEND"          "EXONSEQSTART"        "GCCONTENT"
# #  [9] "GENEBIOTYPE"         "GENEID"              "GENEIDVERSION"       "GENENAME"
# # [13] "GENESEQEND"          "GENESEQSTART"        "INTERPROACCESSION"   "ISCIRCULAR"
# # [17] "PROTDOMEND"          "PROTDOMSTART"        "PROTEINDOMAINID"     "PROTEINDOMAINSOURCE"
# # [21] "PROTEINID"           "PROTEINSEQUENCE"     "SEQCOORDSYSTEM"      "SEQLENGTH"
# # [25] "SEQNAME"             "SEQSTRAND"           "SYMBOL"              "TXBIOTYPE"
# # [29] "TXCDSSEQEND"         "TXCDSSEQSTART"       "TXEXTERNALNAME"      "TXID"
# # [33] "TXIDVERSION"         "TXISCANONICAL"       "TXNAME"              "TXSEQEND"
# # [37] "TXSEQSTART"          "TXSUPPORTLEVEL"      "UNIPROTDB"           "UNIPROTID"
# # [41] "UNIPROTMAPPINGTYPE"
# > keytypes(db)
# #  [1] "ENTREZID"            "EXONID"              "GENEBIOTYPE"         "GENEID"
# #  [5] "GENENAME"            "PROTDOMID"           "PROTEINDOMAINID"     "PROTEINDOMAINSOURCE"
# #  [9] "PROTEINID"           "SEQNAME"             "SEQSTRAND"           "SYMBOL"
# # [13] "TXBIOTYPE"           "TXID"                "TXNAME"              "UNIPROTID"
#
#
# select(db, keys=ids, columns=c("UNIPROTID", "GENEID","SYMBOL"), keytype="UNIPROTID")
#
#
#
# library(EnsDb.Hsapiens.v86)
# db <- EnsDb.Hsapiens.v86
# tbl.idMap <- select(db, keys=ids, columns=c("UNIPROTID", "GENEID","SYMBOL"), keytype="UNIPROTID")
# map <- mapIds(EnsDb.Hsapiens.v79, rownames(tbl), "SYMBOL", "GENEID")
#
#
