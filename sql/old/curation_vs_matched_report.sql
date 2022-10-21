select Modification.id,
       category,
       descriptor,
       el_name_short,
       el_name_long,
       modification_type,
       notes,
       size_bp,
       group_concat(distinct Mcgs.curated_gene_symbols)     as curated_gene_symbols,
       curated_protein_name,
       curated_enzyme_name,
       group_concat(distinct PSmgse.match_gene_symbols_etc) as matched_gene_symbols_etc,
       group_concat(distinct PSmn.match_names)              as matched_names,
       group_concat(distinct PSua.uniprot_accessions)       as matched_uniprot_accessions,
       group_concat(distinct PSoa.other_accessions)         as matched_other_accessions,
       group_concat(distinct PSgtl.go_term_labels)          as matched_go_term_labels
from Modification
         left join Modification_curated_gene_symbols Mcgs on Modification.id = Mcgs.Modification_id
         left join PartsSequence PS on Modification.id = PS.associated_part
         left join PartsSequence_go_term_labels PSgtl on PS.id = PSgtl.PartsSequence_id
         left join PartsSequence_other_accessions PSoa on PS.id = PSoa.PartsSequence_id
         left join PartsSequence_uniprot_accessions PSua on PS.id = PSua.PartsSequence_id
         left join PartsSequence_match_ec_numbers PSmen on PS.id = PSmen.PartsSequence_id
         left join PartsSequence_match_gene_symbols_etc PSmgse on PS.id = PSmgse.PartsSequence_id
         left join PartsSequence_match_names PSmn on PS.id = PSmn.PartsSequence_id
         left join PartsSequence_scientific_names PSsn on PS.id = PSsn.PartsSequence_id
-- left join PartsS
where Mcgs.curated_gene_symbols is not null
group by Modification.id,
         category,
         descriptor,
         el_name_short,
         el_name_long,
         modification_type,
         notes,
         size_bp,
         curated_protein_name,
         curated_enzyme_name;