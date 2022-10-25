

CREATE TABLE "Database" (
	modification_set TEXT, 
	organism_set TEXT, 
	person_set TEXT, 
	sequence_set TEXT, 
	strain_set TEXT, 
	PRIMARY KEY (modification_set, organism_set, person_set, sequence_set, strain_set)
);

CREATE TABLE "Organism" (
	comment TEXT, 
	id TEXT NOT NULL, 
	name TEXT, 
	species_name TEXT, 
	strain_agnostic_taxid TEXT, 
	strain_value TEXT, 
	abbreviation TEXT, 
	PRIMARY KEY (id)
);

CREATE TABLE "Person" (
	id TEXT NOT NULL, 
	date_joined DATETIME, 
	email TEXT, 
	first_name TEXT, 
	is_staff BOOLEAN NOT NULL, 
	is_superuser BOOLEAN NOT NULL, 
	last_name TEXT, 
	username TEXT NOT NULL, 
	PRIMARY KEY (id)
);

CREATE TABLE "Modification" (
	id TEXT NOT NULL, 
	aa_change TEXT, 
	bio_safety_level VARCHAR(7) NOT NULL, 
	category VARCHAR(21), 
	creator TEXT NOT NULL, 
	descriptor VARCHAR(22), 
	el_name_long TEXT NOT NULL, 
	el_name_short TEXT, 
	element_organism TEXT, 
	modification_type VARCHAR(26), 
	modifications_genes TEXT, 
	notes TEXT, 
	position TEXT, 
	principal_investigator TEXT NOT NULL, 
	size_bp INTEGER, 
	status VARCHAR(11) NOT NULL, 
	subcategory_size TEXT, 
	curated_protein_name TEXT, 
	curated_enzyme_name TEXT, 
	curated_uniprot_accession TEXT, 
	PRIMARY KEY (id), 
	FOREIGN KEY(creator) REFERENCES "Person" (id), 
	FOREIGN KEY(principal_investigator) REFERENCES "Person" (id)
);

CREATE TABLE "Strain" (
	bio_safety_level VARCHAR(7), 
	creator TEXT, 
	funding_source VARCHAR(11), 
	genotype_phenotype TEXT, 
	host_organism TEXT, 
	id TEXT NOT NULL, 
	intellectual_property TEXT, 
	keywords TEXT, 
	name TEXT, 
	notes TEXT, 
	principal_investigator TEXT, 
	"references" TEXT, 
	status VARCHAR(11), 
	summary TEXT, 
	PRIMARY KEY (id), 
	FOREIGN KEY(creator) REFERENCES "Person" (id), 
	FOREIGN KEY(host_organism) REFERENCES "Organism" (id), 
	FOREIGN KEY(principal_investigator) REFERENCES "Person" (id)
);

CREATE TABLE "PartsSequence" (
	id TEXT NOT NULL, 
	associated_part TEXT, 
	date_added TEXT, 
	nt_sequence TEXT, 
	seq_name TEXT, 
	seq_type VARCHAR(9), 
	PRIMARY KEY (id), 
	FOREIGN KEY(associated_part) REFERENCES "Modification" (id)
);

CREATE TABLE "Modification_curated_gene_symbols" (
	backref_id TEXT, 
	curated_gene_symbols TEXT, 
	PRIMARY KEY (backref_id, curated_gene_symbols), 
	FOREIGN KEY(backref_id) REFERENCES "Modification" (id)
);

CREATE TABLE "Modification_part_ofs" (
	backref_id TEXT, 
	part_ofs TEXT, 
	PRIMARY KEY (backref_id, part_ofs), 
	FOREIGN KEY(backref_id) REFERENCES "Modification" (id)
);

CREATE TABLE "Strain_biosample_accessions" (
	backref_id TEXT, 
	biosample_accessions TEXT, 
	PRIMARY KEY (backref_id, biosample_accessions), 
	FOREIGN KEY(backref_id) REFERENCES "Strain" (id)
);

CREATE TABLE "Strain_external_urls" (
	backref_id TEXT, 
	external_urls TEXT, 
	PRIMARY KEY (backref_id, external_urls), 
	FOREIGN KEY(backref_id) REFERENCES "Strain" (id)
);

CREATE TABLE "Strain_genome_accessions" (
	backref_id TEXT, 
	genome_accessions TEXT, 
	PRIMARY KEY (backref_id, genome_accessions), 
	FOREIGN KEY(backref_id) REFERENCES "Strain" (id)
);

CREATE TABLE "Strain_has_parts" (
	backref_id TEXT, 
	has_parts TEXT, 
	PRIMARY KEY (backref_id, has_parts), 
	FOREIGN KEY(backref_id) REFERENCES "Strain" (id)
);

CREATE TABLE "Strain_selection_markers" (
	backref_id TEXT, 
	selection_markers TEXT, 
	PRIMARY KEY (backref_id, selection_markers), 
	FOREIGN KEY(backref_id) REFERENCES "Strain" (id)
);

CREATE TABLE "PartsSequence_go_term_ids" (
	backref_id TEXT, 
	go_term_ids TEXT, 
	PRIMARY KEY (backref_id, go_term_ids), 
	FOREIGN KEY(backref_id) REFERENCES "PartsSequence" (id)
);

CREATE TABLE "PartsSequence_go_term_labels" (
	backref_id TEXT, 
	go_term_labels TEXT, 
	PRIMARY KEY (backref_id, go_term_labels), 
	FOREIGN KEY(backref_id) REFERENCES "PartsSequence" (id)
);

CREATE TABLE "PartsSequence_match_ec_numbers" (
	backref_id TEXT, 
	match_ec_numbers TEXT, 
	PRIMARY KEY (backref_id, match_ec_numbers), 
	FOREIGN KEY(backref_id) REFERENCES "PartsSequence" (id)
);

CREATE TABLE "PartsSequence_match_gene_symbols_etc" (
	backref_id TEXT, 
	match_gene_symbols_etc TEXT, 
	PRIMARY KEY (backref_id, match_gene_symbols_etc), 
	FOREIGN KEY(backref_id) REFERENCES "PartsSequence" (id)
);

CREATE TABLE "PartsSequence_match_names" (
	backref_id TEXT, 
	match_names TEXT, 
	PRIMARY KEY (backref_id, match_names), 
	FOREIGN KEY(backref_id) REFERENCES "PartsSequence" (id)
);

CREATE TABLE "PartsSequence_other_accessions" (
	backref_id TEXT, 
	other_accessions TEXT, 
	PRIMARY KEY (backref_id, other_accessions), 
	FOREIGN KEY(backref_id) REFERENCES "PartsSequence" (id)
);

CREATE TABLE "PartsSequence_scientific_names" (
	backref_id TEXT, 
	scientific_names TEXT, 
	PRIMARY KEY (backref_id, scientific_names), 
	FOREIGN KEY(backref_id) REFERENCES "PartsSequence" (id)
);

CREATE TABLE "PartsSequence_taxon_ids" (
	backref_id TEXT, 
	taxon_ids TEXT, 
	PRIMARY KEY (backref_id, taxon_ids), 
	FOREIGN KEY(backref_id) REFERENCES "PartsSequence" (id)
);

CREATE TABLE "PartsSequence_uniprot_accessions" (
	backref_id TEXT, 
	uniprot_accessions TEXT, 
	PRIMARY KEY (backref_id, uniprot_accessions), 
	FOREIGN KEY(backref_id) REFERENCES "PartsSequence" (id)
);
