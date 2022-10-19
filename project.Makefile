# no annotation for non-coding sequences yet
# uniprot annotations for GFPs not especially good
# doesn't yet do any special handling of deletion flanks
# may require ssh tunel to postgres

#max_eval=1e-20
max_eval=1e-5
blast_thread_count=11
#blast_thread_count=15

RUN=poetry run
SCHEMA_DIR = src/synbio_schema/schema/
SCHEMA_NAME = synbio_schema
SCHEMA_EXTENSION = .yaml
SCHEMA_FILE = $(SCHEMA_DIR)$(SCHEMA_NAME)$(SCHEMA_EXTENSION)

.PHONY: \
blast_res_to_sqlite \
confirm_invalid \
interval_clustering \
jsonschema_validation \
project_all \
project_clean \
squeaky_clean \
uniprot_approach

project_all: project_clean jsonschema_validation resources/felix_dump.db \
target/seq2ids.fasta blastdbs/swissprot.psq taxdb.bti target/seq2ids_v_uniprot.tsv \
data/fpbase.fasta data/fpbase.fasta.psq target/seq2ids_v_fpbase.tsv \
blast_res_to_sqlite resources/swiss_entries.json interval_clustering \
resources/synbio_database.json resources/synbio_database.db

squeaky_clean: project_clean
	rm -rf blastdbs/*
	rm -rf data/*
	rm -rf resources/swiss_entries.yaml
	rm -rf swissprot*
	rm -rf taxdb*

project_clean:
	rm -rf resources/*.db
	rm -rf resources/*.json
	rm -rf resources/*.yaml
	rm -rf resources/felix_dump.db
	rm -rf resources/linting_log.tsv
	rm -rf target/*

resources/%.json: $(SCHEMA_FILE) src/data/examples/%.yaml
	$(RUN) linkml-convert \
		--output $@ \
		--target-class Person \
		--schema $^

resources/%.schema.json: $(SCHEMA_FILE)
	$(RUN) gen-json-schema $< > $@

jsonschema_validation: resources/person.json resources/$(SCHEMA_NAME).schema.json
	$(RUN) jsonschema -i $^

#resources/person_set_database.json: $(SCHEMA_FILE) src/data/examples/person_set_database.yaml
#	$(RUN) linkml-convert \
#		--output $@ \
#		--target-class Database \
#		--schema $^
#	$(RUN) jsonschema -i $@ resources/synbio_schema.schema.json
#
#resources/invalid_person_set_database.json: $(SCHEMA_FILE) src/data/examples/invalid_person_set_database.yaml
#	$(RUN) linkml-convert \
#		--output $@ \
#		--target-class Database \
#		--schema $^
#	$(RUN) jsonschema -i $@ resources/synbio_schema.schema.json

#jv2:resources/person_set_database.json resources/$(SCHEMA_NAME).schema.json
#	$(RUN) jsonschema -i $^
#
## can't run this as is because the yaml can't be converted. provide invalid json instead
#confirm_invalid: src/data/examples/invalid_person.json resources/$(SCHEMA_NAME).schema.json
#	! $(RUN) jsonschema -i $^

# --- the following rules will only work for users who have access to the felix database on either the bicoid or staufer servers ---

# skipping tables that don't have at least two populated rows with two populated columns:
#   the entities that one would expect to find ion some of these tables can actually be found in other tables, like...
#
#   arabidopsis_seeds, parts_comments, proteins, strains
# also skipping biological_samples and it's linked tables for this iteration
# sqlite3.OperationalError: table sqlite_master may not be modified messages do not seem to indicate that a table's rows were not copied to SQLite
# but foreign keys do seem to get lost
# todo add more connection parameterization into utils/pgsql2sqlite.sh
# $1 = dbuser
# $2 = port
# $3 = comma-separated list of tables to migrate
# $4 = output file's base name
# omitted: password
# still hardcoded: hostname, dbname, directory and extension for output

resources/felix_dump.db:
	- $(RUN) bash utils/pgsql2sqlite.sh mam 1111 \
	parts,auth_user,species,parts_accessions,parts_sequences,modifications,modifications_genes,selection_markers,plasmids,parts_parameters,external_urls,sub_parts \
	$(basename $(notdir $@))
	sqlite3 $@ "update auth_user set password = NULL;"

# SQLITE -> FASTA
target/seq2ids.fasta: resources/felix_dump.db
	mkdir -p target
	$(RUN) seq2ids \
		--min_len 51 \
		--max_len 50000 \
		--sqlite_file $< \
		--fasta_out $@ \
		--metadata_tsv_out $(subst .fasta,.tsv,$@)

# download NCBI's pre-indexed swissprot BLAST database
# probably doesn't include trembl, so decreased search space
blastdbs/swissprot.psq:
	mkdir -p blastdbs
	wget https://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz
	tar -xzvf swissprot.tar.gz --directory blastdbs

# get NCBI's taxon ID supplement for blast searches
# subsequent steps seem to expect it at the root of the repo
taxdb.bti: blastdbs/swissprot.psq
	wget https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
	tar -xzvf taxdb.tar.gz --directory .

# blast the FELIX proteins against swissprot.
#  did we limit this to insertions somewhere? or some other constraints?
# output is technically misnamed because the search database was swissprot not uniprot
# assumes blastx is on the path
# https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
target/seq2ids_v_uniprot.tsv: target/seq2ids.fasta taxdb.bti
	blastx \
		-query $< \
		-db blastdbs/swissprot \
		-num_threads ${blast_thread_count} \
		-out target/at_delim_blast.txt \
		-evalue ${max_eval} \
		-outfmt "6 delim=@ qacc qcovhsp qcovs qstart qend bitscore score evalue length pident sacc sstart send sallacc sseqid sallseqid stitle salltitles staxids sskingdoms sscinames sblastnames scomnames"
	cat target/at_delim_blast.txt | tr '@' '\t' > target/tab_delim_blast.tsv
	$(RUN) add_col \
		--tsv_in target/tab_delim_blast.tsv \
		--tsv_out $@ \
		--col_val swissprot
	rm target/at_delim_blast.txt target/tab_delim_blast.tsv

# get fluorescent protein database sequences from API
# accessing additional annotations might require GraphQL query in addition to REST API
data/fpbase.fasta:
	mkdir -p data
	$(RUN) get_fluor_prot_seqs \
		--fasta_out data/fpbase.fasta \
		--json_out resources/fpbase.json \
		--yaml_out resources/fpbase.yaml

# index fluorescent protein database for blast
data/fpbase.fasta.psq: data/fpbase.fasta
	makeblastdb -in $< -dbtype prot

target/seq2ids_v_fpbase.tsv: target/seq2ids.fasta data/fpbase.fasta.psq taxdb.bti resources/felix_dump.db
	blastx \
		-query $< \
		-db data/fpbase.fasta \
		-num_threads ${blast_thread_count} \
		-out target/at_delim_blast.txt \
		-evalue ${max_eval} \
		-outfmt "6 delim=@ qacc qcovhsp qcovs qstart qend bitscore score evalue length pident sacc sstart send sallacc sseqid sallseqid stitle salltitles staxids sskingdoms sscinames sblastnames scomnames"
	cat target/at_delim_blast.txt | tr '@' '\t' > target/tab_delim_blast.tsv
	$(RUN) add_col \
		--tsv_in target/tab_delim_blast.tsv \
		--tsv_out $@ \
		--col_val fpbase
	rm target/at_delim_blast.txt target/tab_delim_blast.tsv

blast_res_to_sqlite: target/seq2ids_v_uniprot.tsv target/seq2ids_v_fpbase.tsv
	sqlite3 resources/felix_dump.db < sql/create_blast_results_table.sql
	sqlite3 resources/felix_dump.db ".mode tabs" ".import target/seq2ids_v_uniprot.tsv blast_results" ""
	sqlite3 resources/felix_dump.db ".mode tabs" ".import target/seq2ids_v_fpbase.tsv blast_results" ""
	sqlite3 resources/felix_dump.db < sql/without_attachement.sql
	sqlite3 resources/felix_dump.db < sql/indices.sql

resources/swiss_entries.json:
# blast_res_to_sqlite
	$(RUN) get_uniprot_entries \
		--sqlite_db_fp resources/felix_dump.db \
		--swiss_entries_dump_fp resources/swiss_entries.json \
		--blast_results_table blast_results \
		--requests_cache_name uniprot_entries_cache

interval_clustering:
	# todo is sleeping really necessary?
	sleep 60
	$(RUN) python src/synbio_schema/scripts/interval_clustering.py

resources/synbio_database.yaml:
	$(RUN) python src/synbio_schema/scripts/celniker2yaml.py

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