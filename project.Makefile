# not retrieving annotations for non-coding sequences yet
# not retrieving annotations for flanking sequences yet
# uniprot annotations for GFPs not especially good, so additionally using fpbase
# doesn't yet do any special handling of deletion flanks
# may require ssh tunel to postgres for fetching LBNL_FELIX_TRACKING_DB contents

#blast_thread_count=15
#max_eval=1e-20
BLAST_THREAD_COUNT=11
MAX_E_VAL=1e-5

LBNL_FELIX_TRACKING_DB = resources/lbnl_felix_tracking_dump.db
# using taxdb to get taxon IDs might be necessary for some databases
#BLAST_OUTFMT = "6 delim=@ qacc qcovhsp qcovs qstart qend bitscore score evalue length pident sacc sstart send sallacc sseqid sallseqid stitle salltitles staxids sskingdoms sscinames sblastnames scomnames"
BLAST_OUTFMT = "6 delim=@ qacc qcovhsp qcovs qstart qend bitscore score evalue length pident sacc sstart send sallacc sseqid sallseqid stitle salltitles"


RUN=poetry run
SCHEMA_DIR = src/synbio_schema/schema/
SCHEMA_EXTENSION = .yaml
SCHEMA_NAME = synbio_schema
SCHEMA_FILE = $(SCHEMA_DIR)$(SCHEMA_NAME)$(SCHEMA_EXTENSION)
# src/synbio_schema/schema/synbio_schema.yaml

.PHONY: \
confirm_invalid \
jsonschema_validation \
project_all \
lbnl_felix_clean \
squeaky_clean

#project_all: \

completion_flags/lbnl_felix_tracking_dump: \
lbnl_felix_clean \
jsonschema_validation \
resources/lbnl_felix_tracking_dump.db \
completion_flags/all_blast_databases \
resources/seq2ids_v_uniprot.tsv \
resources/seq2ids_v_fpbase.tsv \
completion_flags/blast_res_to_sqlite \
resources/swiss_entries.json \
completion_flags/interval_clustering \
resources/synbio_database.yaml

squeaky_clean: lbnl_felix_clean clean_blast_databases
	echo $(SCHEMA_FILE)
	rm -rf completion_flags/*
	mkdir -p completion_flags

clean_blast_databases:
	rm -rf downloads/swissprot.tar.gz
	rm -rf blastdbs/*
	rm -rf completion_flags/all_blast_databases
	mkdir -p blastdbs

lbnl_felix_clean:
	rm -rf resources/*.db
	rm -rf resources/*.fasta
	rm -rf resources/*.json
	rm -rf resources/*.tsv
	rm -rf resources/*.yaml
	mkdir -p resources

completion_flags/all_blast_databases: \
blastdbs/swissprot.psq \
blastdbs/fpbase.fasta.psq
	echo "done" > $@

resources/%.json: src/synbio_schema/schema/synbio_schema.yaml src/data/examples/%.yaml
	$(RUN) linkml-convert \
		--output $@ \
		--target-class Person \
		--schema $^

resources/%.schema.json: src/synbio_schema/schema/synbio_schema.yaml
	$(RUN) gen-json-schema $< > $@

# this just a basic validation of the schema against known good data
jsonschema_validation: resources/person.json resources/$(SCHEMA_NAME).schema.json
	$(RUN) jsonschema -i $^

# --- the following rules will only work for users who have access to the felix database on either the bicoid or staufer servers ---

# skipping tables that don't have at least two populated rows with two populated columns:
#   the entities that one would expect to find ion some of these tables can actually be found in other tables, like...
#   arabidopsis_seeds, parts_comments, proteins, strains
# also skipping biological_samples and its linked tables for this iteration
# sqlite3.OperationalError: table sqlite_master may not be modified messages do not seem to indicate that a table's rows were not copied to SQLite
# but foreign keys do seem to get lost
# todo add more connection parameterization into utils/pgsql2sqlite.sh
# $1 = dbuser
# $2 = port
# $3 = comma-separated list of tables to migrate
# $4 = output file's base name
# omitted: password
# still hardcoded: hostname, dbname, directory and extension for output

resources/lbnl_felix_tracking_dump.db:
	- $(RUN) bash utils/pg2sqlite.sh mam 1111 \
	parts,auth_user,species,parts_accessions,parts_sequences,modifications,modifications_genes,selection_markers,plasmids,parts_parameters,external_urls,sub_parts \
	$(basename $(notdir $@))
	sqlite3 $@ "update auth_user set password = NULL;"

# SQLITE -> FASTA
# 		--max_len 50000
resources/seq2ids.fasta: resources/lbnl_felix_tracking_dump.db
	mkdir -p target
	$(RUN) seq2ids \
		--min_len 51 \
		--max_len 500 \
		--sqlite_file $< \
		--fasta_out $@ \
		--metadata_tsv_out $(subst .fasta,.tsv,$@)

# download NCBI's pre-indexed swissprot BLAST database
# doesn't include trembl, so the search space is limited
# but it does complete quickly
downloads/swissprot.tar.gz:
	mkdir -p blastdbs
	mkdir -p downloads
	wget https://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz --directory-prefix $(dir $@)
	#tar -xzvf swissprot.tar.gz --directory blastdbs

blastdbs/swissprot.psq: downloads/swissprot.tar.gz
	tar -xzvf $< --directory blastdbs

## get NCBI's taxon ID supplement for blast searches
## subsequent steps seem to expect it at the root of the repo
## todo does this add anything compared to the taxon info from the swissprot entity lookup?
#taxdb.bti: blastdbs/swissprot.psq
#	wget https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
#	tar -xzvf taxdb.tar.gz --directory .

# blast the FELIX proteins against swissprot.
#  make sure we didn't limit this to modifications or sequences of type "insertion"s somewhere? or some other constraints?
# output is technically misnamed because the search database was swissprot not uniprot
# assumes blastx is on the path
# https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
# 		-out target/at_delim_blast.txt \

resources/seq2ids_v_uniprot.tsv: resources/seq2ids.fasta blastdbs/swissprot.psq
	blastx \
		-query $< \
		-db $(basename $(word 2, $^)) \
		-num_threads $(BLAST_THREAD_COUNT) \
		-evalue $(MAX_E_VAL) \
		-outfmt $(BLAST_OUTFMT) | tr '@' '\t' > $@.temp
	$(RUN) add_col \
		--tsv_in $@.temp \
		--tsv_out $@ \
		--col_val swissprot
	rm $@.temp

# get fluorescent protein database sequences from API
# accessing additional annotations might require GraphQL query in addition to REST API
blastdbs/fpbase.%:
	mkdir -p data
	$(RUN) get_fluor_prot_seqs \
		--fasta_out $(basename $@).fasta \
		--json_out $(basename $@).json \
		--yaml_out $(basename $@).yaml

# index fluorescent protein database for blast
blastdbs/fpbase.fasta.psq: blastdbs/fpbase.fasta
	makeblastdb -in $< -dbtype prot

resources/seq2ids_v_fpbase.tsv: resources/seq2ids.fasta blastdbs/fpbase.fasta.psq
	blastx \
		-query $< \
		-db $(basename $(word 2, $^)) \
		-num_threads $(BLAST_THREAD_COUNT) \
		-evalue $(MAX_E_VAL) \
		-outfmt $(BLAST_OUTFMT) | tr '@' '\t' > $@.temp
	$(RUN) add_col \
		--tsv_in $@.temp \
		--tsv_out $@ \
		--col_val fpbase
	rm $@.temp

# resources/seq2ids_v_fpbase.tsv restarts
# resources/seq2ids_v_uniprot.tsv restarts
completion_flags/blast_res_to_sqlite: resources/lbnl_felix_tracking_dump.db
	sqlite3 $< < sql/create_blast_results_table.sql
	sqlite3 $<  ".mode tabs" ".import resources/seq2ids_v_uniprot.tsv blast_results" ""
	sqlite3 $<  ".mode tabs" ".import resources/seq2ids_v_fpbase.tsv blast_results" ""
	sqlite3 $<  < sql/calc_seq_lens.sql
	sqlite3 $<  < sql/indices.sql
	echo "done" >  $@


# completion_flags/blast_res_to_sqlite restarts
resources/swiss_entries.json: resources/lbnl_felix_tracking_dump.db completion_flags/blast_res_to_sqlite
	$(RUN) get_uniprot_entries \
		--blast_results_table blast_results \
		--requests_cache_name uniprot_entries_cache \
		--sqlite_db_fp $< \
		--swiss_entries_dump_fp $@

# sqlite3.OperationalError: no such table: uniprot_annotation_scores
completion_flags/interval_clustering: \
completion_flags/blast_res_to_sqlite \
resources/lbnl_felix_tracking_dump.db \
resources/swiss_entries.json
	$(RUN) interval_clustering \
		--sqlite_file $(word 2, $^)
	echo "done" > $@

#

# completion_flags/interval_clustering restarts
# resources/swiss_entries.json restarts
resources/synbio_database.yaml: blastdbs/fpbase.json resources/lbnl_felix_tracking_dump.db resources/swiss_entries.json completion_flags/interval_clustering
	$(RUN) lbnl_felix_tracking2linkml \
		--fpbase_entries_fp $(word 1, $^) \
		--sqlite_input_fp $(word 2, $^) \
		--swiss_entries_fp $(word 3, $^) \
		--yaml_out $@

resources/synbio_database.json: src/synbio_schema/schema/synbio_schema.yaml resources/synbio_database.yaml
	$(RUN) linkml-convert \
		--output $@ \
		--target-class Database \
		--schema $^

resources/synbio_database.db: src/synbio_schema/schema/synbio_schema.yaml resources/synbio_database.json
	$(RUN) linkml-sqldb dump \
		--db $@ \
		--target-class Database \
		--schema $^

## can't run this as is because the yaml can't be converted. provide invalid json instead
#confirm_invalid: src/data/examples/invalid_person.json resources/$(SCHEMA_NAME).schema.json
#	! $(RUN) jsonschema -i $^