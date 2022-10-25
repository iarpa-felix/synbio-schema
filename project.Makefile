# not retrieving annotations for non-coding sequences yet
# not retrieving annotations for flanking sequences yet
# doesn't yet do any special handling of deletion flanks
# uniprot annotations for GFPs not especially good, so additionally using fpbase
# may require ssh tunel to postgres for fetching LBNL_FELIX_TRACKING_DB contents

BLAST_THREAD_COUNT=11
MAX_E_VAL=1e-5

LBNL_FELIX_TRACKING_DB = resources/lbnl_felix_tracking_dump.db
BLAST_OUTFMT = "6 delim=@ qacc qcovhsp qcovs qstart qend bitscore score evalue length pident sacc sstart send sallacc sseqid sallseqid stitle salltitles"

RUN=poetry run
SCHEMA_DIR = src/synbio_schema/schema/
SCHEMA_EXTENSION = .yaml
SCHEMA_NAME = synbio_schema
SCHEMA_FILE = $(SCHEMA_DIR)$(SCHEMA_NAME)$(SCHEMA_EXTENSION)

.PHONY: \
confirm_invalid \
jsonschema_validation \
lbnl_felix_clean

completion_flags/lbnl_felix_tracking_dump: \
lbnl_felix_clean \
jsonschema_validation \
completion_flags/all_blast_databases \
resources/seq2ids_v_uniprot.tsv \
resources/seq2ids_v_fpbase.tsv \
completion_flags/blast_res_to_sqlite \
resources/swiss_entries.json \
completion_flags/interval_clustering \
resources/synbio_database.yaml \
resources/synbio_database.db \
resources/synbio_database.json \
resources/synbio_database.ttl

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
	rm -rf resources/*.ttl
	rm -rf resources/*.yaml
	mkdir -p resources
	rm -rf completion_flags/*
	mkdir -p completion_flags

completion_flags/all_blast_databases: \
blastdbs/swissprot.psq \
blastdbs/fpbase.fasta.psq
	touch $@

# ---

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

# ---

# --- the following rules will only work for users who have access to the felix database on either the bicoid or staufer servers ---

# sqlite3.OperationalError: table sqlite_master may not be modified messages do not seem to indicate that a table's rows were not copied to SQLite
# but foreign keys do seem to get lost
# todo add more connection parameterization into utils/pgsql2sqlite.sh
# $1 = dbuser
# $2 = port
# $3 = comma-separated list of tables to migrate
# $4 = output file's base name
# omitted: password
# still hardcoded: hostname, dbname, directory and extension for output

completion_flags/pg2sqlite:
	- $(RUN) bash utils/pg2sqlite.sh mam 1111 \
	parts,auth_user,species,parts_accessions,parts_sequences,modifications,modifications_genes,selection_markers,plasmids,parts_parameters,external_urls,sub_parts \
	lbnl_felix_tracking_dump
	sqlite3 resources/lbnl_felix_tracking_dump.db "update auth_user set password = NULL;"
	touch $@

# SQLITE -> FASTA
# 		--max_len 50000
resources/seq2ids.fasta: completion_flags/pg2sqlite
	$(RUN) seq2ids \
		--min_len 51 \
		--max_len 500 \
		--sqlite_file resources/lbnl_felix_tracking_dump.db \
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
	touch $@


# blast the FELIX proteins against swissprot.
#  make sure we didn't limit this to modifications or sequences of type "insertion"s somewhere? or some other constraints?
# output is technically misnamed because the search database was swissprot not uniprot
# assumes blastx is on the path
# https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
# 		-out target/at_delim_blast.txt \

resources/seq2ids_v_uniprot.tsv: \
resources/seq2ids.fasta \
blastdbs/swissprot.psq
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

completion_flags/blast_res_to_sqlite:  completion_flags/pg2sqlite resources/seq2ids_v_fpbase.tsv resources/seq2ids_v_uniprot.tsv
	sqlite3 resources/lbnl_felix_tracking_dump.db < sql/create_blast_results_table.sql
	sqlite3 resources/lbnl_felix_tracking_dump.db  ".mode tabs" ".import resources/seq2ids_v_uniprot.tsv blast_results" ""
	sqlite3 resources/lbnl_felix_tracking_dump.db  ".mode tabs" ".import resources/seq2ids_v_fpbase.tsv blast_results" ""
	sqlite3 resources/lbnl_felix_tracking_dump.db  < sql/calc_seq_lens.sql
	sqlite3 resources/lbnl_felix_tracking_dump.db  < sql/indices.sql
	touch  $@


resources/swiss_entries.json: \
completion_flags/blast_res_to_sqlite
	$(RUN) get_uniprot_entries \
		--blast_results_table blast_results \
		--requests_cache_name uniprot_entries_cache \
		--sqlite_db_fp resources/lbnl_felix_tracking_dump.db \
		--swiss_entries_dump_fp $@


completion_flags/interval_clustering: \
completion_flags/blast_res_to_sqlite \
resources/swiss_entries.json
	$(RUN) interval_clustering \
		--sqlite_file resources/lbnl_felix_tracking_dump.db
	echo "done" > $@


resources/synbio_database.yaml: \
completion_flags/interval_clustering \
blastdbs/fpbase.json \
resources/lbnl_felix_tracking_dump.db \
resources/swiss_entries.json
	$(RUN) lbnl_felix_tracking2linkml \
		--fpbase_entries_fp $(word 2, $^) \
		--sqlite_input_fp $(word 3, $^) \
		--swiss_entries_fp $(word 4, $^) \
		--yaml_out $@

resources/synbio_database.json: src/synbio_schema/schema/synbio_schema.yaml resources/synbio_database.yaml
	$(RUN) linkml-convert \
		--output $@ \
		--target-class Database \
		--schema $^

resources/synbio_database.db: src/synbio_schema/schema/synbio_schema.yaml resources/synbio_database.yaml
	$(RUN) linkml-sqldb dump \
		--db $@ \
		--target-class Database \
		--schema $^

## can't run this as is because the yaml can't be converted. provide invalid json instead
#confirm_invalid: src/data/examples/invalid_person.json resources/$(SCHEMA_NAME).schema.json
#	! $(RUN) jsonschema -i $^

resources/raw_species.tsv: resources/lbnl_felix_tracking_dump.db
	sqlite3 -readonly -header -csv -separator '	' $< "select * from species;" > resources/raw_species.tsv

#resources/raw_species.yaml: resources/raw_species.tsv
#	$(RUN) schemauto generalize-csv $< -o $@

resources/synbio_database.ttl: src/synbio_schema/schema/synbio_schema.yaml resources/synbio_database.yaml
	$(RUN) linkml-convert \
		--output $@ \
		--target-class Database \
		--schema $^