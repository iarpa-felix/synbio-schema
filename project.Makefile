# no annotation for non-coding sequences yet
# uniprot annotations for GFPs not especially good

RUN=poetry run
SCHEMA_DIR = src/synbio_schema/schema/
SCHEMA_NAME = synbio_schema
SCHEMA_EXTENSION = .yaml
SCHEMA_FILE = $(SCHEMA_DIR)$(SCHEMA_NAME)$(SCHEMA_EXTENSION)

.PHONY: project_clean project_all jsonschema_validation confirm_invalid

# probably don't want to release this with resources/felix_dump.db as part of project_all
project_all: project_clean  jsonschema_validation resources/felix_dump.db

# resources/linting_log.tsv

project_clean:
	rm -rf resources/linting_log.tsv
	rm -rf resources/*.db
	rm -rf resources/*.yaml
	rm -rf resources/*.json

#  -c, --config FILE
#  -f, --format [terminal|markdown|json|tsv]
#  -o, --output FILENAME
#  --fix / --no-fix

## or make lint
#resources/linting_log.tsv: $(SCHEMA_FILE)
#	$(RUN) linkml-lint \
#		--config resources/config/linter_config_default.yaml \
#		--format tsv  \
#		--output $@ $<

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
	- $(RUN) sh utils/pgsql2sqlite.sh mam 1111 \
	parts,auth_user,species,parts_accessions,parts_sequences,modifications,selection_markers,plasmids,parts_parameters,external_urls,sub_parts \
	$(basename $(notdir $@))
	sqlite3 $@ "update auth_user set password = NULL;"

