# Auto generated from synbio_schema.yaml by pythongen.py version: 0.9.0
# Generation date: 2022-10-25T12:09:24
# Schema: synbio
#
# id: http://www.semanticweb.org/mam/ontologies/2022/7/synbio
# description: synbio
# license: https://creativecommons.org/publicdomain/zero/1.0/

import dataclasses
import sys
import re
from jsonasobj2 import JsonObj, as_dict
from typing import Optional, List, Union, Dict, ClassVar, Any
from dataclasses import dataclass
from linkml_runtime.linkml_model.meta import EnumDefinition, PermissibleValue, PvFormulaOptions

from linkml_runtime.utils.slot import Slot
from linkml_runtime.utils.metamodelcore import empty_list, empty_dict, bnode
from linkml_runtime.utils.yamlutils import YAMLRoot, extended_str, extended_float, extended_int
from linkml_runtime.utils.dataclass_extensions_376 import dataclasses_init_fn_with_kwargs
from linkml_runtime.utils.formatutils import camelcase, underscore, sfx
from linkml_runtime.utils.enumerations import EnumDefinitionImpl
from rdflib import Namespace, URIRef
from linkml_runtime.utils.curienamespace import CurieNamespace
from linkml_runtime.linkml_model.types import Boolean, Datetime, Integer, String, Uriorcurie
from linkml_runtime.utils.metamodelcore import Bool, URIorCURIE, XSDDateTime

metamodel_version = "1.7.0"
version = None

# Overwrite dataclasses _init_fn to add **kwargs in __init__
dataclasses._init_fn = dataclasses_init_fn_with_kwargs

# Namespaces
GO = CurieNamespace('GO', 'http://amigo.geneontology.org/amigo/term/GO:')
IAO = CurieNamespace('IAO', 'http://purl.obolibrary.org/obo/IAO_')
IF = CurieNamespace('IF', 'https://example.com/modification/')
MS = CurieNamespace('MS', 'https://example.com/strain/')
NCBITAXON = CurieNamespace('NCBITaxon', 'http://purl.obolibrary.org/obo/NCBITaxon_')
NCIT = CurieNamespace('NCIT', 'http://purl.obolibrary.org/obo/NCIT_')
OBI = CurieNamespace('OBI', 'http://purl.obolibrary.org/obo/OBI_')
SO = CurieNamespace('SO', 'http://purl.obolibrary.org/obo/SO_')
DCTERMS = CurieNamespace('dcterms', 'http://purl.org/dc/terms/')
LINKML = CurieNamespace('linkml', 'https://w3id.org/linkml/')
ORGANISM = CurieNamespace('organism', 'https://example.com/organism/')
PERSON = CurieNamespace('person', 'https://example.com/person/')
SEQUENCE = CurieNamespace('sequence', 'https://example.com/sequence/')
SYNBIO = CurieNamespace('synbio', 'https://example.com/synbio/')
XSD = CurieNamespace('xsd', 'http://www.w3.org/2001/XMLSchema#')
DEFAULT_ = SYNBIO


# Types

# Class references
class PartsSequenceId(extended_str):
    pass


class StrainId(extended_str):
    pass


class ModificationId(extended_str):
    pass


class OrganismId(extended_str):
    pass


class PersonId(extended_str):
    pass


@dataclass
class PartsSequence(YAMLRoot):
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = SYNBIO.PartsSequence
    class_class_curie: ClassVar[str] = "synbio:PartsSequence"
    class_name: ClassVar[str] = "PartsSequence"
    class_model_uri: ClassVar[URIRef] = SYNBIO.PartsSequence

    id: Union[str, PartsSequenceId] = None
    associated_part: Optional[Union[str, ModificationId]] = None
    date_added: Optional[str] = None
    go_term_ids: Optional[Union[Union[str, URIorCURIE], List[Union[str, URIorCURIE]]]] = empty_list()
    go_term_labels: Optional[Union[str, List[str]]] = empty_list()
    match_ec_numbers: Optional[Union[str, List[str]]] = empty_list()
    match_gene_symbols_etc: Optional[Union[str, List[str]]] = empty_list()
    match_names: Optional[Union[str, List[str]]] = empty_list()
    nt_sequence: Optional[str] = None
    other_accessions: Optional[Union[str, List[str]]] = empty_list()
    scientific_names: Optional[Union[str, List[str]]] = empty_list()
    seq_name: Optional[str] = None
    seq_type: Optional[Union[str, "SeqTypeEnum"]] = None
    taxon_ids: Optional[Union[str, List[str]]] = empty_list()
    uniprot_accessions: Optional[Union[str, List[str]]] = empty_list()

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, PartsSequenceId):
            self.id = PartsSequenceId(self.id)

        if self.associated_part is not None and not isinstance(self.associated_part, ModificationId):
            self.associated_part = ModificationId(self.associated_part)

        if self.date_added is not None and not isinstance(self.date_added, str):
            self.date_added = str(self.date_added)

        if not isinstance(self.go_term_ids, list):
            self.go_term_ids = [self.go_term_ids] if self.go_term_ids is not None else []
        self.go_term_ids = [v if isinstance(v, URIorCURIE) else URIorCURIE(v) for v in self.go_term_ids]

        if not isinstance(self.go_term_labels, list):
            self.go_term_labels = [self.go_term_labels] if self.go_term_labels is not None else []
        self.go_term_labels = [v if isinstance(v, str) else str(v) for v in self.go_term_labels]

        if not isinstance(self.match_ec_numbers, list):
            self.match_ec_numbers = [self.match_ec_numbers] if self.match_ec_numbers is not None else []
        self.match_ec_numbers = [v if isinstance(v, str) else str(v) for v in self.match_ec_numbers]

        if not isinstance(self.match_gene_symbols_etc, list):
            self.match_gene_symbols_etc = [self.match_gene_symbols_etc] if self.match_gene_symbols_etc is not None else []
        self.match_gene_symbols_etc = [v if isinstance(v, str) else str(v) for v in self.match_gene_symbols_etc]

        if not isinstance(self.match_names, list):
            self.match_names = [self.match_names] if self.match_names is not None else []
        self.match_names = [v if isinstance(v, str) else str(v) for v in self.match_names]

        if self.nt_sequence is not None and not isinstance(self.nt_sequence, str):
            self.nt_sequence = str(self.nt_sequence)

        if not isinstance(self.other_accessions, list):
            self.other_accessions = [self.other_accessions] if self.other_accessions is not None else []
        self.other_accessions = [v if isinstance(v, str) else str(v) for v in self.other_accessions]

        if not isinstance(self.scientific_names, list):
            self.scientific_names = [self.scientific_names] if self.scientific_names is not None else []
        self.scientific_names = [v if isinstance(v, str) else str(v) for v in self.scientific_names]

        if self.seq_name is not None and not isinstance(self.seq_name, str):
            self.seq_name = str(self.seq_name)

        if self.seq_type is not None and not isinstance(self.seq_type, SeqTypeEnum):
            self.seq_type = SeqTypeEnum(self.seq_type)

        if not isinstance(self.taxon_ids, list):
            self.taxon_ids = [self.taxon_ids] if self.taxon_ids is not None else []
        self.taxon_ids = [v if isinstance(v, str) else str(v) for v in self.taxon_ids]

        if not isinstance(self.uniprot_accessions, list):
            self.uniprot_accessions = [self.uniprot_accessions] if self.uniprot_accessions is not None else []
        self.uniprot_accessions = [v if isinstance(v, str) else str(v) for v in self.uniprot_accessions]

        super().__post_init__(**kwargs)


@dataclass
class Strain(YAMLRoot):
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = SYNBIO.Strain
    class_class_curie: ClassVar[str] = "synbio:Strain"
    class_name: ClassVar[str] = "Strain"
    class_model_uri: ClassVar[URIRef] = SYNBIO.Strain

    id: Union[str, StrainId] = None
    bio_safety_level: Optional[Union[str, "BioSafetyLevelEnum"]] = None
    biosample_accessions: Optional[Union[str, List[str]]] = empty_list()
    creator: Optional[Union[str, PersonId]] = None
    external_urls: Optional[Union[str, List[str]]] = empty_list()
    funding_source: Optional[Union[str, "FundingSourceEnum"]] = None
    genome_accessions: Optional[Union[str, List[str]]] = empty_list()
    genotype_phenotype: Optional[str] = None
    has_parts: Optional[Union[str, List[str]]] = empty_list()
    host_organism: Optional[Union[str, OrganismId]] = None
    intellectual_property: Optional[str] = None
    keywords: Optional[str] = None
    name: Optional[str] = None
    notes: Optional[str] = None
    principal_investigator: Optional[Union[str, PersonId]] = None
    references: Optional[str] = None
    selection_markers: Optional[Union[str, List[str]]] = empty_list()
    status: Optional[Union[str, "StatusEnum"]] = None
    summary: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, StrainId):
            self.id = StrainId(self.id)

        if self.bio_safety_level is not None and not isinstance(self.bio_safety_level, BioSafetyLevelEnum):
            self.bio_safety_level = BioSafetyLevelEnum(self.bio_safety_level)

        if not isinstance(self.biosample_accessions, list):
            self.biosample_accessions = [self.biosample_accessions] if self.biosample_accessions is not None else []
        self.biosample_accessions = [v if isinstance(v, str) else str(v) for v in self.biosample_accessions]

        if self.creator is not None and not isinstance(self.creator, PersonId):
            self.creator = PersonId(self.creator)

        if not isinstance(self.external_urls, list):
            self.external_urls = [self.external_urls] if self.external_urls is not None else []
        self.external_urls = [v if isinstance(v, str) else str(v) for v in self.external_urls]

        if self.funding_source is not None and not isinstance(self.funding_source, FundingSourceEnum):
            self.funding_source = FundingSourceEnum(self.funding_source)

        if not isinstance(self.genome_accessions, list):
            self.genome_accessions = [self.genome_accessions] if self.genome_accessions is not None else []
        self.genome_accessions = [v if isinstance(v, str) else str(v) for v in self.genome_accessions]

        if self.genotype_phenotype is not None and not isinstance(self.genotype_phenotype, str):
            self.genotype_phenotype = str(self.genotype_phenotype)

        if not isinstance(self.has_parts, list):
            self.has_parts = [self.has_parts] if self.has_parts is not None else []
        self.has_parts = [v if isinstance(v, str) else str(v) for v in self.has_parts]

        if self.host_organism is not None and not isinstance(self.host_organism, OrganismId):
            self.host_organism = OrganismId(self.host_organism)

        if self.intellectual_property is not None and not isinstance(self.intellectual_property, str):
            self.intellectual_property = str(self.intellectual_property)

        if self.keywords is not None and not isinstance(self.keywords, str):
            self.keywords = str(self.keywords)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.notes is not None and not isinstance(self.notes, str):
            self.notes = str(self.notes)

        if self.principal_investigator is not None and not isinstance(self.principal_investigator, PersonId):
            self.principal_investigator = PersonId(self.principal_investigator)

        if self.references is not None and not isinstance(self.references, str):
            self.references = str(self.references)

        if not isinstance(self.selection_markers, list):
            self.selection_markers = [self.selection_markers] if self.selection_markers is not None else []
        self.selection_markers = [v if isinstance(v, str) else str(v) for v in self.selection_markers]

        if self.status is not None and not isinstance(self.status, StatusEnum):
            self.status = StatusEnum(self.status)

        if self.summary is not None and not isinstance(self.summary, str):
            self.summary = str(self.summary)

        super().__post_init__(**kwargs)


@dataclass
class Database(YAMLRoot):
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = SYNBIO.Database
    class_class_curie: ClassVar[str] = "synbio:Database"
    class_name: ClassVar[str] = "Database"
    class_model_uri: ClassVar[URIRef] = SYNBIO.Database

    modification_set: Optional[Union[Dict[Union[str, ModificationId], Union[dict, "Modification"]], List[Union[dict, "Modification"]]]] = empty_dict()
    organism_set: Optional[Union[Dict[Union[str, OrganismId], Union[dict, "Organism"]], List[Union[dict, "Organism"]]]] = empty_dict()
    person_set: Optional[Union[Dict[Union[str, PersonId], Union[dict, "Person"]], List[Union[dict, "Person"]]]] = empty_dict()
    sequence_set: Optional[Union[Dict[Union[str, PartsSequenceId], Union[dict, PartsSequence]], List[Union[dict, PartsSequence]]]] = empty_dict()
    strain_set: Optional[Union[Dict[Union[str, StrainId], Union[dict, Strain]], List[Union[dict, Strain]]]] = empty_dict()

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        self._normalize_inlined_as_list(slot_name="modification_set", slot_type=Modification, key_name="id", keyed=True)

        self._normalize_inlined_as_list(slot_name="organism_set", slot_type=Organism, key_name="id", keyed=True)

        self._normalize_inlined_as_list(slot_name="person_set", slot_type=Person, key_name="id", keyed=True)

        self._normalize_inlined_as_list(slot_name="sequence_set", slot_type=PartsSequence, key_name="id", keyed=True)

        self._normalize_inlined_as_list(slot_name="strain_set", slot_type=Strain, key_name="id", keyed=True)

        super().__post_init__(**kwargs)


@dataclass
class Modification(YAMLRoot):
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = SYNBIO.Modification
    class_class_curie: ClassVar[str] = "synbio:Modification"
    class_name: ClassVar[str] = "Modification"
    class_model_uri: ClassVar[URIRef] = SYNBIO.Modification

    id: Union[str, ModificationId] = None
    bio_safety_level: Union[str, "BioSafetyLevelEnum"] = None
    creator: Union[str, PersonId] = None
    el_name_long: str = None
    principal_investigator: Union[str, PersonId] = None
    status: Union[str, "StatusEnum"] = None
    aa_change: Optional[str] = None
    category: Optional[Union[str, "CategoryEnum"]] = None
    descriptor: Optional[Union[str, "DescriptorEnum"]] = None
    el_name_short: Optional[str] = None
    element_organism: Optional[str] = None
    modification_type: Optional[Union[str, "ModificationTypeEnum"]] = None
    modifications_genes: Optional[str] = None
    notes: Optional[str] = None
    position: Optional[str] = None
    size_bp: Optional[int] = None
    subcategory_size: Optional[str] = None
    curated_gene_symbols: Optional[Union[str, List[str]]] = empty_list()
    curated_protein_name: Optional[str] = None
    curated_enzyme_name: Optional[str] = None
    curated_uniprot_accession: Optional[str] = None
    part_ofs: Optional[Union[str, List[str]]] = empty_list()

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, ModificationId):
            self.id = ModificationId(self.id)

        if self._is_empty(self.bio_safety_level):
            self.MissingRequiredField("bio_safety_level")
        if not isinstance(self.bio_safety_level, BioSafetyLevelEnum):
            self.bio_safety_level = BioSafetyLevelEnum(self.bio_safety_level)

        if self._is_empty(self.creator):
            self.MissingRequiredField("creator")
        if not isinstance(self.creator, PersonId):
            self.creator = PersonId(self.creator)

        if self._is_empty(self.el_name_long):
            self.MissingRequiredField("el_name_long")
        if not isinstance(self.el_name_long, str):
            self.el_name_long = str(self.el_name_long)

        if self._is_empty(self.principal_investigator):
            self.MissingRequiredField("principal_investigator")
        if not isinstance(self.principal_investigator, PersonId):
            self.principal_investigator = PersonId(self.principal_investigator)

        if self._is_empty(self.status):
            self.MissingRequiredField("status")
        if not isinstance(self.status, StatusEnum):
            self.status = StatusEnum(self.status)

        if self.aa_change is not None and not isinstance(self.aa_change, str):
            self.aa_change = str(self.aa_change)

        if self.category is not None and not isinstance(self.category, CategoryEnum):
            self.category = CategoryEnum(self.category)

        if self.descriptor is not None and not isinstance(self.descriptor, DescriptorEnum):
            self.descriptor = DescriptorEnum(self.descriptor)

        if self.el_name_short is not None and not isinstance(self.el_name_short, str):
            self.el_name_short = str(self.el_name_short)

        if self.element_organism is not None and not isinstance(self.element_organism, str):
            self.element_organism = str(self.element_organism)

        if self.modification_type is not None and not isinstance(self.modification_type, ModificationTypeEnum):
            self.modification_type = ModificationTypeEnum(self.modification_type)

        if self.modifications_genes is not None and not isinstance(self.modifications_genes, str):
            self.modifications_genes = str(self.modifications_genes)

        if self.notes is not None and not isinstance(self.notes, str):
            self.notes = str(self.notes)

        if self.position is not None and not isinstance(self.position, str):
            self.position = str(self.position)

        if self.size_bp is not None and not isinstance(self.size_bp, int):
            self.size_bp = int(self.size_bp)

        if self.subcategory_size is not None and not isinstance(self.subcategory_size, str):
            self.subcategory_size = str(self.subcategory_size)

        if not isinstance(self.curated_gene_symbols, list):
            self.curated_gene_symbols = [self.curated_gene_symbols] if self.curated_gene_symbols is not None else []
        self.curated_gene_symbols = [v if isinstance(v, str) else str(v) for v in self.curated_gene_symbols]

        if self.curated_protein_name is not None and not isinstance(self.curated_protein_name, str):
            self.curated_protein_name = str(self.curated_protein_name)

        if self.curated_enzyme_name is not None and not isinstance(self.curated_enzyme_name, str):
            self.curated_enzyme_name = str(self.curated_enzyme_name)

        if self.curated_uniprot_accession is not None and not isinstance(self.curated_uniprot_accession, str):
            self.curated_uniprot_accession = str(self.curated_uniprot_accession)

        if not isinstance(self.part_ofs, list):
            self.part_ofs = [self.part_ofs] if self.part_ofs is not None else []
        self.part_ofs = [v if isinstance(v, str) else str(v) for v in self.part_ofs]

        super().__post_init__(**kwargs)


@dataclass
class Organism(YAMLRoot):
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = SYNBIO.Organism
    class_class_curie: ClassVar[str] = "synbio:Organism"
    class_name: ClassVar[str] = "Organism"
    class_model_uri: ClassVar[URIRef] = SYNBIO.Organism

    id: Union[str, OrganismId] = None
    comment: Optional[str] = None
    name: Optional[str] = None
    species_name: Optional[str] = None
    strain_agnostic_taxid: Optional[str] = None
    strain_value: Optional[str] = None
    abbreviation: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, OrganismId):
            self.id = OrganismId(self.id)

        if self.comment is not None and not isinstance(self.comment, str):
            self.comment = str(self.comment)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.species_name is not None and not isinstance(self.species_name, str):
            self.species_name = str(self.species_name)

        if self.strain_agnostic_taxid is not None and not isinstance(self.strain_agnostic_taxid, str):
            self.strain_agnostic_taxid = str(self.strain_agnostic_taxid)

        if self.strain_value is not None and not isinstance(self.strain_value, str):
            self.strain_value = str(self.strain_value)

        if self.abbreviation is not None and not isinstance(self.abbreviation, str):
            self.abbreviation = str(self.abbreviation)

        super().__post_init__(**kwargs)


@dataclass
class Person(YAMLRoot):
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = SYNBIO.Person
    class_class_curie: ClassVar[str] = "synbio:Person"
    class_name: ClassVar[str] = "Person"
    class_model_uri: ClassVar[URIRef] = SYNBIO.Person

    id: Union[str, PersonId] = None
    is_staff: Union[bool, Bool] = None
    is_superuser: Union[bool, Bool] = None
    username: str = None
    date_joined: Optional[Union[str, XSDDateTime]] = None
    email: Optional[str] = None
    first_name: Optional[str] = None
    last_name: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, PersonId):
            self.id = PersonId(self.id)

        if self._is_empty(self.is_staff):
            self.MissingRequiredField("is_staff")
        if not isinstance(self.is_staff, Bool):
            self.is_staff = Bool(self.is_staff)

        if self._is_empty(self.is_superuser):
            self.MissingRequiredField("is_superuser")
        if not isinstance(self.is_superuser, Bool):
            self.is_superuser = Bool(self.is_superuser)

        if self._is_empty(self.username):
            self.MissingRequiredField("username")
        if not isinstance(self.username, str):
            self.username = str(self.username)

        if self.date_joined is not None and not isinstance(self.date_joined, XSDDateTime):
            self.date_joined = XSDDateTime(self.date_joined)

        if self.email is not None and not isinstance(self.email, str):
            self.email = str(self.email)

        if self.first_name is not None and not isinstance(self.first_name, str):
            self.first_name = str(self.first_name)

        if self.last_name is not None and not isinstance(self.last_name, str):
            self.last_name = str(self.last_name)

        super().__post_init__(**kwargs)


# Enumerations
class SeqTypeEnum(EnumDefinitionImpl):

    deletion = PermissibleValue(text="deletion")
    sequence = PermissibleValue(text="sequence")
    insertion = PermissibleValue(text="insertion")
    flank1 = PermissibleValue(text="flank1")
    flank2 = PermissibleValue(text="flank2")

    _defn = EnumDefinition(
        name="SeqTypeEnum",
    )

class IntellectualPropertyEnum(EnumDefinitionImpl):

    _defn = EnumDefinition(
        name="IntellectualPropertyEnum",
    )

    @classmethod
    def _addvals(cls):
        setattr(cls, "Limited License Agreement required from Scarab Genomics",
                PermissibleValue(text="Limited License Agreement required from Scarab Genomics") )

class FundingSourceEnum(EnumDefinitionImpl):

    _defn = EnumDefinition(
        name="FundingSourceEnum",
    )

    @classmethod
    def _addvals(cls):
        setattr(cls, "IARPA-FELIX",
                PermissibleValue(text="IARPA-FELIX") )

class StatusEnum(EnumDefinitionImpl):

    Abandoned = PermissibleValue(text="Abandoned")

    _defn = EnumDefinition(
        name="StatusEnum",
    )

    @classmethod
    def _addvals(cls):
        setattr(cls, "In Progress",
                PermissibleValue(text="In Progress") )

class ModificationTypeEnum(EnumDefinitionImpl):

    deletion = PermissibleValue(text="deletion")
    frameshift = PermissibleValue(text="frameshift")
    insertion = PermissibleValue(text="insertion")
    inversion = PermissibleValue(text="inversion")
    plasmid = PermissibleValue(text="plasmid")
    substitution = PermissibleValue(text="substitution")
    transition = PermissibleValue(text="transition")
    CDSpartial = PermissibleValue(text="CDSpartial")
    reassortment = PermissibleValue(text="reassortment")
    subtitutions = PermissibleValue(text="subtitutions")

    _defn = EnumDefinition(
        name="ModificationTypeEnum",
    )

    @classmethod
    def _addvals(cls):
        setattr(cls, "amber stop codon",
                PermissibleValue(text="amber stop codon") )
        setattr(cls, "insertion site",
                PermissibleValue(text="insertion site") )
        setattr(cls, "ochre stop codon",
                PermissibleValue(text="ochre stop codon") )
        setattr(cls, "203TAGTAGGTACT",
                PermissibleValue(text="203TAGTAGGTACT") )
        setattr(cls, "492ACTTAAGCCAGAAAATTTA",
                PermissibleValue(text="492ACTTAAGCCAGAAAATTTA") )
        setattr(cls, "T474A-A490T-T492A",
                PermissibleValue(text="T474A-A490T-T492A") )
        setattr(cls, "chrXIII-chrII",
                PermissibleValue(text="chrXIII-chrII") )
        setattr(cls, "compound deletion",
                PermissibleValue(text="compound deletion") )
        setattr(cls, "compound insertion",
                PermissibleValue(text="compound insertion") )
        setattr(cls, "fluorescent/epitope",
                PermissibleValue(text="fluorescent/epitope") )
        setattr(cls, "in-frame deletion",
                PermissibleValue(text="in-frame deletion") )
        setattr(cls, "plasmid element",
                PermissibleValue(text="plasmid element") )
        setattr(cls, "substitution transition",
                PermissibleValue(text="substitution transition") )
        setattr(cls, "substitution transversion",
                PermissibleValue(text="substitution transversion") )
        setattr(cls, "transition plasmid element",
                PermissibleValue(text="transition plasmid element") )

class BioSafetyLevelEnum(EnumDefinitionImpl):

    _defn = EnumDefinition(
        name="BioSafetyLevelEnum",
    )

    @classmethod
    def _addvals(cls):
        setattr(cls, "Level 1",
                PermissibleValue(text="Level 1") )
        setattr(cls, "Level 2",
                PermissibleValue(text="Level 2") )

class CategoryEnum(EnumDefinitionImpl):

    biosynthetic = PermissibleValue(text="biosynthetic")
    compound = PermissibleValue(text="compound")
    GOI = PermissibleValue(text="GOI")
    promoter = PermissibleValue(text="promoter")
    terminator = PermissibleValue(text="terminator")

    _defn = EnumDefinition(
        name="CategoryEnum",
    )

    @classmethod
    def _addvals(cls):
        setattr(cls, "antibiotic resistance",
                PermissibleValue(text="antibiotic resistance") )
        setattr(cls, "control element",
                PermissibleValue(text="control element") )
        setattr(cls, "fluorescent/epitope",
                PermissibleValue(text="fluorescent/epitope") )
        setattr(cls, "origin/structural",
                PermissibleValue(text="origin/structural") )
        setattr(cls, "toxin/virulence",
                PermissibleValue(text="toxin/virulence") )

class DescriptorEnum(EnumDefinitionImpl):

    CDS = PermissibleValue(text="CDS")
    CDSpartial = PermissibleValue(text="CDSpartial")
    eGFP = PermissibleValue(text="eGFP")
    enhancer = PermissibleValue(text="enhancer")
    epitope = PermissibleValue(text="epitope")
    fragment = PermissibleValue(text="fragment")
    G = PermissibleValue(text="G")
    gene = PermissibleValue(text="gene")
    opt = PermissibleValue(text="opt")
    promoter = PermissibleValue(text="promoter")
    scar = PermissibleValue(text="scar")
    SUMO = PermissibleValue(text="SUMO")
    terminator = PermissibleValue(text="terminator")
    VP16 = PermissibleValue(text="VP16")
    GFP = PermissibleValue(text="GFP")

    _defn = EnumDefinition(
        name="DescriptorEnum",
    )

    @classmethod
    def _addvals(cls):
        setattr(cls, "3pUTR",
                PermissibleValue(text="3pUTR") )
        setattr(cls, "5pUTR",
                PermissibleValue(text="5pUTR") )
        setattr(cls, "activiation domain",
                PermissibleValue(text="activiation domain") )
        setattr(cls, "CDSpartial-3prime",
                PermissibleValue(text="CDSpartial-3prime") )
        setattr(cls, "CDSpartial-5prime",
                PermissibleValue(text="CDSpartial-5prime") )
        setattr(cls, "insertion site",
                PermissibleValue(text="insertion site") )
        setattr(cls, "multiple cloning site",
                PermissibleValue(text="multiple cloning site") )
        setattr(cls, "ribosomal binding site",
                PermissibleValue(text="ribosomal binding site") )
        setattr(cls, "cleavage site",
                PermissibleValue(text="cleavage site") )
        setattr(cls, "self-cleaving-peptide",
                PermissibleValue(text="self-cleaving-peptide") )

# Slots
class slots:
    pass

slots.genome_accessions = Slot(uri=SYNBIO.genome_accessions, name="genome_accessions", curie=SYNBIO.curie('genome_accessions'),
                   model_uri=SYNBIO.genome_accessions, domain=None, range=Optional[Union[str, List[str]]])

slots.biosample_accessions = Slot(uri=SYNBIO.biosample_accessions, name="biosample_accessions", curie=SYNBIO.curie('biosample_accessions'),
                   model_uri=SYNBIO.biosample_accessions, domain=None, range=Optional[Union[str, List[str]]])

slots.selection_markers = Slot(uri=SYNBIO.selection_markers, name="selection_markers", curie=SYNBIO.curie('selection_markers'),
                   model_uri=SYNBIO.selection_markers, domain=None, range=Optional[Union[str, List[str]]])

slots.external_urls = Slot(uri=SYNBIO.external_urls, name="external_urls", curie=SYNBIO.curie('external_urls'),
                   model_uri=SYNBIO.external_urls, domain=None, range=Optional[Union[str, List[str]]])

slots.curated_uniprot_accession = Slot(uri=SYNBIO.curated_uniprot_accession, name="curated_uniprot_accession", curie=SYNBIO.curie('curated_uniprot_accession'),
                   model_uri=SYNBIO.curated_uniprot_accession, domain=None, range=Optional[str])

slots.has_parts = Slot(uri=SYNBIO.has_parts, name="has_parts", curie=SYNBIO.curie('has_parts'),
                   model_uri=SYNBIO.has_parts, domain=None, range=Optional[Union[str, List[str]]])

slots.part_ofs = Slot(uri=SYNBIO.part_ofs, name="part_ofs", curie=SYNBIO.curie('part_ofs'),
                   model_uri=SYNBIO.part_ofs, domain=None, range=Optional[Union[str, List[str]]])

slots.curated_protein_name = Slot(uri=SYNBIO.curated_protein_name, name="curated_protein_name", curie=SYNBIO.curie('curated_protein_name'),
                   model_uri=SYNBIO.curated_protein_name, domain=None, range=Optional[str])

slots.curated_enzyme_name = Slot(uri=SYNBIO.curated_enzyme_name, name="curated_enzyme_name", curie=SYNBIO.curie('curated_enzyme_name'),
                   model_uri=SYNBIO.curated_enzyme_name, domain=None, range=Optional[str])

slots.curated_gene_symbols = Slot(uri=SYNBIO.curated_gene_symbols, name="curated_gene_symbols", curie=SYNBIO.curie('curated_gene_symbols'),
                   model_uri=SYNBIO.curated_gene_symbols, domain=None, range=Optional[str])

slots.other_accessions = Slot(uri=SYNBIO.other_accessions, name="other_accessions", curie=SYNBIO.curie('other_accessions'),
                   model_uri=SYNBIO.other_accessions, domain=None, range=Optional[Union[str, List[str]]])

slots.match_gene_symbols_etc = Slot(uri=SYNBIO.match_gene_symbols_etc, name="match_gene_symbols_etc", curie=SYNBIO.curie('match_gene_symbols_etc'),
                   model_uri=SYNBIO.match_gene_symbols_etc, domain=None, range=Optional[Union[str, List[str]]])

slots.match_ec_numbers = Slot(uri=SYNBIO.match_ec_numbers, name="match_ec_numbers", curie=SYNBIO.curie('match_ec_numbers'),
                   model_uri=SYNBIO.match_ec_numbers, domain=None, range=Optional[Union[str, List[str]]])

slots.match_names = Slot(uri=SYNBIO.match_names, name="match_names", curie=SYNBIO.curie('match_names'),
                   model_uri=SYNBIO.match_names, domain=None, range=Optional[Union[str, List[str]]])

slots.scientific_names = Slot(uri=SYNBIO.scientific_names, name="scientific_names", curie=SYNBIO.curie('scientific_names'),
                   model_uri=SYNBIO.scientific_names, domain=None, range=Optional[Union[str, List[str]]])

slots.taxon_ids = Slot(uri=SYNBIO.taxon_ids, name="taxon_ids", curie=SYNBIO.curie('taxon_ids'),
                   model_uri=SYNBIO.taxon_ids, domain=None, range=Optional[Union[str, List[str]]])

slots.uniprot_accessions = Slot(uri=SYNBIO.uniprot_accessions, name="uniprot_accessions", curie=SYNBIO.curie('uniprot_accessions'),
                   model_uri=SYNBIO.uniprot_accessions, domain=None, range=Optional[Union[str, List[str]]])

slots.go_term_ids = Slot(uri=SYNBIO.go_term_ids, name="go_term_ids", curie=SYNBIO.curie('go_term_ids'),
                   model_uri=SYNBIO.go_term_ids, domain=None, range=Optional[Union[Union[str, URIorCURIE], List[Union[str, URIorCURIE]]]])

slots.go_term_labels = Slot(uri=SYNBIO.go_term_labels, name="go_term_labels", curie=SYNBIO.curie('go_term_labels'),
                   model_uri=SYNBIO.go_term_labels, domain=None, range=Optional[Union[str, List[str]]])

slots.associated_part = Slot(uri=SYNBIO.associated_part, name="associated_part", curie=SYNBIO.curie('associated_part'),
                   model_uri=SYNBIO.associated_part, domain=None, range=Optional[Union[str, ModificationId]])

slots.seq_name = Slot(uri=SYNBIO.seq_name, name="seq_name", curie=SYNBIO.curie('seq_name'),
                   model_uri=SYNBIO.seq_name, domain=None, range=Optional[str])

slots.seq_type = Slot(uri=SYNBIO.seq_type, name="seq_type", curie=SYNBIO.curie('seq_type'),
                   model_uri=SYNBIO.seq_type, domain=None, range=Optional[Union[str, "SeqTypeEnum"]])

slots.date_added = Slot(uri=SYNBIO.date_added, name="date_added", curie=SYNBIO.curie('date_added'),
                   model_uri=SYNBIO.date_added, domain=None, range=Optional[str])

slots.nt_sequence = Slot(uri=SYNBIO.nt_sequence, name="nt_sequence", curie=SYNBIO.curie('nt_sequence'),
                   model_uri=SYNBIO.nt_sequence, domain=None, range=Optional[str])

slots.host_organism = Slot(uri=SYNBIO.host_organism, name="host_organism", curie=SYNBIO.curie('host_organism'),
                   model_uri=SYNBIO.host_organism, domain=None, range=Optional[Union[str, OrganismId]])

slots.strain_set = Slot(uri=SYNBIO.strain_set, name="strain_set", curie=SYNBIO.curie('strain_set'),
                   model_uri=SYNBIO.strain_set, domain=None, range=Optional[Union[Dict[Union[str, StrainId], Union[dict, Strain]], List[Union[dict, Strain]]]])

slots.sequence_set = Slot(uri=SYNBIO.sequence_set, name="sequence_set", curie=SYNBIO.curie('sequence_set'),
                   model_uri=SYNBIO.sequence_set, domain=None, range=Optional[Union[Dict[Union[str, PartsSequenceId], Union[dict, PartsSequence]], List[Union[dict, PartsSequence]]]])

slots.name = Slot(uri=SYNBIO.name, name="name", curie=SYNBIO.curie('name'),
                   model_uri=SYNBIO.name, domain=None, range=Optional[str])

slots.genotype_phenotype = Slot(uri=SYNBIO.genotype_phenotype, name="genotype_phenotype", curie=SYNBIO.curie('genotype_phenotype'),
                   model_uri=SYNBIO.genotype_phenotype, domain=None, range=Optional[str])

slots.intellectual_property = Slot(uri=SYNBIO.intellectual_property, name="intellectual_property", curie=SYNBIO.curie('intellectual_property'),
                   model_uri=SYNBIO.intellectual_property, domain=None, range=Optional[str])

slots.keywords = Slot(uri=SYNBIO.keywords, name="keywords", curie=SYNBIO.curie('keywords'),
                   model_uri=SYNBIO.keywords, domain=None, range=Optional[str])

slots.summary = Slot(uri=SYNBIO.summary, name="summary", curie=SYNBIO.curie('summary'),
                   model_uri=SYNBIO.summary, domain=None, range=Optional[str])

slots.references = Slot(uri=SYNBIO.references, name="references", curie=SYNBIO.curie('references'),
                   model_uri=SYNBIO.references, domain=None, range=Optional[str])

slots.funding_source = Slot(uri=SYNBIO.funding_source, name="funding_source", curie=SYNBIO.curie('funding_source'),
                   model_uri=SYNBIO.funding_source, domain=None, range=Optional[Union[str, "FundingSourceEnum"]])

slots.modification_set = Slot(uri=SYNBIO.modification_set, name="modification_set", curie=SYNBIO.curie('modification_set'),
                   model_uri=SYNBIO.modification_set, domain=None, range=Optional[Union[Dict[Union[str, ModificationId], Union[dict, Modification]], List[Union[dict, Modification]]]])

slots.category = Slot(uri=SYNBIO.category, name="category", curie=SYNBIO.curie('category'),
                   model_uri=SYNBIO.category, domain=None, range=Optional[Union[str, "CategoryEnum"]])

slots.creator = Slot(uri=SYNBIO.creator, name="creator", curie=SYNBIO.curie('creator'),
                   model_uri=SYNBIO.creator, domain=None, range=Optional[Union[str, PersonId]])

slots.descriptor = Slot(uri=SYNBIO.descriptor, name="descriptor", curie=SYNBIO.curie('descriptor'),
                   model_uri=SYNBIO.descriptor, domain=None, range=Optional[Union[str, "DescriptorEnum"]])

slots.modification_taxon = Slot(uri=SYNBIO.modification_taxon, name="modification_taxon", curie=SYNBIO.curie('modification_taxon'),
                   model_uri=SYNBIO.modification_taxon, domain=None, range=Optional[str])

slots.modification_type = Slot(uri=SYNBIO.modification_type, name="modification_type", curie=SYNBIO.curie('modification_type'),
                   model_uri=SYNBIO.modification_type, domain=None, range=Optional[Union[str, "ModificationTypeEnum"]])

slots.organism_basis = Slot(uri=SYNBIO.organism_basis, name="organism_basis", curie=SYNBIO.curie('organism_basis'),
                   model_uri=SYNBIO.organism_basis, domain=None, range=Optional[str])

slots.organism_set = Slot(uri=SYNBIO.organism_set, name="organism_set", curie=SYNBIO.curie('organism_set'),
                   model_uri=SYNBIO.organism_set, domain=None, range=Optional[Union[Dict[Union[str, OrganismId], Union[dict, Organism]], List[Union[dict, Organism]]]])

slots.person_set = Slot(uri=SYNBIO.person_set, name="person_set", curie=SYNBIO.curie('person_set'),
                   model_uri=SYNBIO.person_set, domain=None, range=Optional[Union[Dict[Union[str, PersonId], Union[dict, Person]], List[Union[dict, Person]]]])

slots.principal_investigator = Slot(uri=SYNBIO.principal_investigator, name="principal_investigator", curie=SYNBIO.curie('principal_investigator'),
                   model_uri=SYNBIO.principal_investigator, domain=None, range=Optional[Union[str, PersonId]])

slots.subcategory_size = Slot(uri=SYNBIO.subcategory_size, name="subcategory_size", curie=SYNBIO.curie('subcategory_size'),
                   model_uri=SYNBIO.subcategory_size, domain=None, range=Optional[str])

slots.license = Slot(uri=DCTERMS.license, name="license", curie=DCTERMS.curie('license'),
                   model_uri=SYNBIO.license, domain=None, range=Optional[str])

slots.aa_change = Slot(uri=SYNBIO.aa_change, name="aa_change", curie=SYNBIO.curie('aa_change'),
                   model_uri=SYNBIO.aa_change, domain=None, range=Optional[str])

slots.abbreviation = Slot(uri=SYNBIO.abbreviation, name="abbreviation", curie=SYNBIO.curie('abbreviation'),
                   model_uri=SYNBIO.abbreviation, domain=None, range=Optional[str])

slots.bio_safety_level = Slot(uri=SYNBIO.bio_safety_level, name="bio_safety_level", curie=SYNBIO.curie('bio_safety_level'),
                   model_uri=SYNBIO.bio_safety_level, domain=None, range=Optional[Union[str, "BioSafetyLevelEnum"]])

slots.comment = Slot(uri=SYNBIO.comment, name="comment", curie=SYNBIO.curie('comment'),
                   model_uri=SYNBIO.comment, domain=None, range=Optional[str])

slots.date_joined = Slot(uri=SYNBIO.date_joined, name="date_joined", curie=SYNBIO.curie('date_joined'),
                   model_uri=SYNBIO.date_joined, domain=None, range=Optional[Union[str, XSDDateTime]])

slots.el_name_long = Slot(uri=SYNBIO.el_name_long, name="el_name_long", curie=SYNBIO.curie('el_name_long'),
                   model_uri=SYNBIO.el_name_long, domain=None, range=Optional[str])

slots.el_name_short = Slot(uri=SYNBIO.el_name_short, name="el_name_short", curie=SYNBIO.curie('el_name_short'),
                   model_uri=SYNBIO.el_name_short, domain=None, range=Optional[str])

slots.element = Slot(uri=SYNBIO.element, name="element", curie=SYNBIO.curie('element'),
                   model_uri=SYNBIO.element, domain=None, range=Optional[str])

slots.element_id = Slot(uri=SYNBIO.element_id, name="element_id", curie=SYNBIO.curie('element_id'),
                   model_uri=SYNBIO.element_id, domain=None, range=Optional[str])

slots.element_name = Slot(uri=SYNBIO.element_name, name="element_name", curie=SYNBIO.curie('element_name'),
                   model_uri=SYNBIO.element_name, domain=None, range=Optional[str])

slots.element_organism = Slot(uri=SYNBIO.element_organism, name="element_organism", curie=SYNBIO.curie('element_organism'),
                   model_uri=SYNBIO.element_organism, domain=None, range=Optional[str])

slots.email = Slot(uri=SYNBIO.email, name="email", curie=SYNBIO.curie('email'),
                   model_uri=SYNBIO.email, domain=None, range=Optional[str])

slots.first_name = Slot(uri=SYNBIO.first_name, name="first_name", curie=SYNBIO.curie('first_name'),
                   model_uri=SYNBIO.first_name, domain=None, range=Optional[str])

slots.id = Slot(uri=SYNBIO.id, name="id", curie=SYNBIO.curie('id'),
                   model_uri=SYNBIO.id, domain=None, range=URIRef)

slots.is_staff = Slot(uri=SYNBIO.is_staff, name="is_staff", curie=SYNBIO.curie('is_staff'),
                   model_uri=SYNBIO.is_staff, domain=None, range=Union[bool, Bool])

slots.is_superuser = Slot(uri=SYNBIO.is_superuser, name="is_superuser", curie=SYNBIO.curie('is_superuser'),
                   model_uri=SYNBIO.is_superuser, domain=None, range=Union[bool, Bool])

slots.last_name = Slot(uri=SYNBIO.last_name, name="last_name", curie=SYNBIO.curie('last_name'),
                   model_uri=SYNBIO.last_name, domain=None, range=Optional[str])

slots.modification_genes = Slot(uri=SYNBIO.modification_genes, name="modification_genes", curie=SYNBIO.curie('modification_genes'),
                   model_uri=SYNBIO.modification_genes, domain=None, range=Optional[str])

slots.modifications_genes = Slot(uri=SYNBIO.modifications_genes, name="modifications_genes", curie=SYNBIO.curie('modifications_genes'),
                   model_uri=SYNBIO.modifications_genes, domain=None, range=Optional[str])

slots.notes = Slot(uri=SYNBIO.notes, name="notes", curie=SYNBIO.curie('notes'),
                   model_uri=SYNBIO.notes, domain=None, range=Optional[str])

slots.position = Slot(uri=SYNBIO.position, name="position", curie=SYNBIO.curie('position'),
                   model_uri=SYNBIO.position, domain=None, range=Optional[str])

slots.size_bp = Slot(uri=SYNBIO.size_bp, name="size_bp", curie=SYNBIO.curie('size_bp'),
                   model_uri=SYNBIO.size_bp, domain=None, range=Optional[int])

slots.species_name = Slot(uri=SYNBIO.species_name, name="species_name", curie=SYNBIO.curie('species_name'),
                   model_uri=SYNBIO.species_name, domain=None, range=Optional[str])

slots.strain_agnostic_taxid = Slot(uri=SYNBIO.strain_agnostic_taxid, name="strain_agnostic_taxid", curie=SYNBIO.curie('strain_agnostic_taxid'),
                   model_uri=SYNBIO.strain_agnostic_taxid, domain=None, range=Optional[str])

slots.status = Slot(uri=SYNBIO.status, name="status", curie=SYNBIO.curie('status'),
                   model_uri=SYNBIO.status, domain=None, range=Optional[Union[str, "StatusEnum"]])

slots.strain_value = Slot(uri=SYNBIO.strain, name="strain_value", curie=SYNBIO.curie('strain'),
                   model_uri=SYNBIO.strain_value, domain=None, range=Optional[str])

slots.username = Slot(uri=SYNBIO.username, name="username", curie=SYNBIO.curie('username'),
                   model_uri=SYNBIO.username, domain=None, range=str)

slots.examples = Slot(uri=LINKML.examples, name="examples", curie=LINKML.curie('examples'),
                   model_uri=SYNBIO.examples, domain=None, range=Optional[str])

slots.generation_date = Slot(uri=LINKML.generation_date, name="generation_date", curie=LINKML.curie('generation_date'),
                   model_uri=SYNBIO.generation_date, domain=None, range=Optional[str])

slots.metamodel_version = Slot(uri=LINKML.metamodel_version, name="metamodel_version", curie=LINKML.curie('metamodel_version'),
                   model_uri=SYNBIO.metamodel_version, domain=None, range=Optional[str])

slots.permissible_values = Slot(uri=LINKML.permissible_values, name="permissible_values", curie=LINKML.curie('permissible_values'),
                   model_uri=SYNBIO.permissible_values, domain=None, range=Optional[str])

slots.source_file = Slot(uri=LINKML.source_file, name="source_file", curie=LINKML.curie('source_file'),
                   model_uri=SYNBIO.source_file, domain=None, range=Optional[str])

slots.source_file_date = Slot(uri=LINKML.source_file_date, name="source_file_date", curie=LINKML.curie('source_file_date'),
                   model_uri=SYNBIO.source_file_date, domain=None, range=Optional[str])

slots.source_file_size = Slot(uri=LINKML.source_file_size, name="source_file_size", curie=LINKML.curie('source_file_size'),
                   model_uri=SYNBIO.source_file_size, domain=None, range=Optional[str])

slots.PartsSequence_id = Slot(uri=SYNBIO.id, name="PartsSequence_id", curie=SYNBIO.curie('id'),
                   model_uri=SYNBIO.PartsSequence_id, domain=PartsSequence, range=Union[str, PartsSequenceId],
                   pattern=re.compile(r'^sequence:\d+$'))

slots.PartsSequence_associated_part = Slot(uri=SYNBIO.associated_part, name="PartsSequence_associated_part", curie=SYNBIO.curie('associated_part'),
                   model_uri=SYNBIO.PartsSequence_associated_part, domain=PartsSequence, range=Optional[Union[str, ModificationId]],
                   pattern=re.compile(r'^IF:\d+$'))

slots.Strain_id = Slot(uri=SYNBIO.id, name="Strain_id", curie=SYNBIO.curie('id'),
                   model_uri=SYNBIO.Strain_id, domain=Strain, range=Union[str, StrainId],
                   pattern=re.compile(r'MS:\d+'))

slots.Database_strain_set = Slot(uri=SYNBIO.strain_set, name="Database_strain_set", curie=SYNBIO.curie('strain_set'),
                   model_uri=SYNBIO.Database_strain_set, domain=Database, range=Optional[Union[Dict[Union[str, StrainId], Union[dict, Strain]], List[Union[dict, Strain]]]])

slots.Database_modification_set = Slot(uri=SYNBIO.modification_set, name="Database_modification_set", curie=SYNBIO.curie('modification_set'),
                   model_uri=SYNBIO.Database_modification_set, domain=Database, range=Optional[Union[Dict[Union[str, ModificationId], Union[dict, "Modification"]], List[Union[dict, "Modification"]]]])

slots.Database_organism_set = Slot(uri=SYNBIO.organism_set, name="Database_organism_set", curie=SYNBIO.curie('organism_set'),
                   model_uri=SYNBIO.Database_organism_set, domain=Database, range=Optional[Union[Dict[Union[str, OrganismId], Union[dict, "Organism"]], List[Union[dict, "Organism"]]]])

slots.Database_person_set = Slot(uri=SYNBIO.person_set, name="Database_person_set", curie=SYNBIO.curie('person_set'),
                   model_uri=SYNBIO.Database_person_set, domain=Database, range=Optional[Union[Dict[Union[str, PersonId], Union[dict, "Person"]], List[Union[dict, "Person"]]]])

slots.Database_sequence_set = Slot(uri=SYNBIO.sequence_set, name="Database_sequence_set", curie=SYNBIO.curie('sequence_set'),
                   model_uri=SYNBIO.Database_sequence_set, domain=Database, range=Optional[Union[Dict[Union[str, PartsSequenceId], Union[dict, PartsSequence]], List[Union[dict, PartsSequence]]]])

slots.Modification_id = Slot(uri=SYNBIO.id, name="Modification_id", curie=SYNBIO.curie('id'),
                   model_uri=SYNBIO.Modification_id, domain=Modification, range=Union[str, ModificationId],
                   pattern=re.compile(r'^IF:\d+$'))

slots.Modification_aa_change = Slot(uri=SYNBIO.aa_change, name="Modification_aa_change", curie=SYNBIO.curie('aa_change'),
                   model_uri=SYNBIO.Modification_aa_change, domain=Modification, range=Optional[str])

slots.Modification_bio_safety_level = Slot(uri=SYNBIO.bio_safety_level, name="Modification_bio_safety_level", curie=SYNBIO.curie('bio_safety_level'),
                   model_uri=SYNBIO.Modification_bio_safety_level, domain=Modification, range=Union[str, "BioSafetyLevelEnum"])

slots.Modification_category = Slot(uri=SYNBIO.category, name="Modification_category", curie=SYNBIO.curie('category'),
                   model_uri=SYNBIO.Modification_category, domain=Modification, range=Optional[Union[str, "CategoryEnum"]])

slots.Modification_creator = Slot(uri=SYNBIO.creator, name="Modification_creator", curie=SYNBIO.curie('creator'),
                   model_uri=SYNBIO.Modification_creator, domain=Modification, range=Union[str, PersonId])

slots.Modification_descriptor = Slot(uri=SYNBIO.descriptor, name="Modification_descriptor", curie=SYNBIO.curie('descriptor'),
                   model_uri=SYNBIO.Modification_descriptor, domain=Modification, range=Optional[Union[str, "DescriptorEnum"]])

slots.Modification_el_name_long = Slot(uri=SYNBIO.el_name_long, name="Modification_el_name_long", curie=SYNBIO.curie('el_name_long'),
                   model_uri=SYNBIO.Modification_el_name_long, domain=Modification, range=str)

slots.Modification_el_name_short = Slot(uri=SYNBIO.el_name_short, name="Modification_el_name_short", curie=SYNBIO.curie('el_name_short'),
                   model_uri=SYNBIO.Modification_el_name_short, domain=Modification, range=Optional[str])

slots.Modification_modification_type = Slot(uri=SYNBIO.modification_type, name="Modification_modification_type", curie=SYNBIO.curie('modification_type'),
                   model_uri=SYNBIO.Modification_modification_type, domain=Modification, range=Optional[Union[str, "ModificationTypeEnum"]])

slots.Modification_notes = Slot(uri=SYNBIO.notes, name="Modification_notes", curie=SYNBIO.curie('notes'),
                   model_uri=SYNBIO.Modification_notes, domain=Modification, range=Optional[str])

slots.Modification_position = Slot(uri=SYNBIO.position, name="Modification_position", curie=SYNBIO.curie('position'),
                   model_uri=SYNBIO.Modification_position, domain=Modification, range=Optional[str])

slots.Modification_principal_investigator = Slot(uri=SYNBIO.principal_investigator, name="Modification_principal_investigator", curie=SYNBIO.curie('principal_investigator'),
                   model_uri=SYNBIO.Modification_principal_investigator, domain=Modification, range=Union[str, PersonId])

slots.Modification_size_bp = Slot(uri=SYNBIO.size_bp, name="Modification_size_bp", curie=SYNBIO.curie('size_bp'),
                   model_uri=SYNBIO.Modification_size_bp, domain=Modification, range=Optional[int])

slots.Modification_status = Slot(uri=SYNBIO.status, name="Modification_status", curie=SYNBIO.curie('status'),
                   model_uri=SYNBIO.Modification_status, domain=Modification, range=Union[str, "StatusEnum"])

slots.Modification_subcategory_size = Slot(uri=SYNBIO.subcategory_size, name="Modification_subcategory_size", curie=SYNBIO.curie('subcategory_size'),
                   model_uri=SYNBIO.Modification_subcategory_size, domain=Modification, range=Optional[str])

slots.Modification_curated_gene_symbols = Slot(uri=SYNBIO.curated_gene_symbols, name="Modification_curated_gene_symbols", curie=SYNBIO.curie('curated_gene_symbols'),
                   model_uri=SYNBIO.Modification_curated_gene_symbols, domain=Modification, range=Optional[Union[str, List[str]]])

slots.Modification_element_organism = Slot(uri=SYNBIO.element_organism, name="Modification_element_organism", curie=SYNBIO.curie('element_organism'),
                   model_uri=SYNBIO.Modification_element_organism, domain=Modification, range=Optional[str])

slots.Modification_modifications_genes = Slot(uri=SYNBIO.modifications_genes, name="Modification_modifications_genes", curie=SYNBIO.curie('modifications_genes'),
                   model_uri=SYNBIO.Modification_modifications_genes, domain=Modification, range=Optional[str])

slots.Organism_id = Slot(uri=SYNBIO.id, name="Organism_id", curie=SYNBIO.curie('id'),
                   model_uri=SYNBIO.Organism_id, domain=Organism, range=Union[str, OrganismId],
                   pattern=re.compile(r'^organism:-?\d+$'))

slots.Organism_abbreviation = Slot(uri=SYNBIO.abbreviation, name="Organism_abbreviation", curie=SYNBIO.curie('abbreviation'),
                   model_uri=SYNBIO.Organism_abbreviation, domain=Organism, range=Optional[str])

slots.Organism_comment = Slot(uri=SYNBIO.comment, name="Organism_comment", curie=SYNBIO.curie('comment'),
                   model_uri=SYNBIO.Organism_comment, domain=Organism, range=Optional[str])

slots.Organism_species_name = Slot(uri=SYNBIO.species_name, name="Organism_species_name", curie=SYNBIO.curie('species_name'),
                   model_uri=SYNBIO.Organism_species_name, domain=Organism, range=Optional[str])

slots.Organism_strain_agnostic_taxid = Slot(uri=SYNBIO.strain_agnostic_taxid, name="Organism_strain_agnostic_taxid", curie=SYNBIO.curie('strain_agnostic_taxid'),
                   model_uri=SYNBIO.Organism_strain_agnostic_taxid, domain=Organism, range=Optional[str])

slots.Organism_strain_value = Slot(uri=SYNBIO.strain, name="Organism_strain_value", curie=SYNBIO.curie('strain'),
                   model_uri=SYNBIO.Organism_strain_value, domain=Organism, range=Optional[str])

slots.Person_id = Slot(uri=SYNBIO.id, name="Person_id", curie=SYNBIO.curie('id'),
                   model_uri=SYNBIO.Person_id, domain=Person, range=Union[str, PersonId],
                   pattern=re.compile(r'^person:\d+$'))

slots.Person_date_joined = Slot(uri=SYNBIO.date_joined, name="Person_date_joined", curie=SYNBIO.curie('date_joined'),
                   model_uri=SYNBIO.Person_date_joined, domain=Person, range=Optional[Union[str, XSDDateTime]])

slots.Person_email = Slot(uri=SYNBIO.email, name="Person_email", curie=SYNBIO.curie('email'),
                   model_uri=SYNBIO.Person_email, domain=Person, range=Optional[str])

slots.Person_first_name = Slot(uri=SYNBIO.first_name, name="Person_first_name", curie=SYNBIO.curie('first_name'),
                   model_uri=SYNBIO.Person_first_name, domain=Person, range=Optional[str])

slots.Person_is_staff = Slot(uri=SYNBIO.is_staff, name="Person_is_staff", curie=SYNBIO.curie('is_staff'),
                   model_uri=SYNBIO.Person_is_staff, domain=Person, range=Union[bool, Bool])

slots.Person_is_superuser = Slot(uri=SYNBIO.is_superuser, name="Person_is_superuser", curie=SYNBIO.curie('is_superuser'),
                   model_uri=SYNBIO.Person_is_superuser, domain=Person, range=Union[bool, Bool])

slots.Person_last_name = Slot(uri=SYNBIO.last_name, name="Person_last_name", curie=SYNBIO.curie('last_name'),
                   model_uri=SYNBIO.Person_last_name, domain=Person, range=Optional[str])

slots.Person_username = Slot(uri=SYNBIO.username, name="Person_username", curie=SYNBIO.curie('username'),
                   model_uri=SYNBIO.Person_username, domain=Person, range=str)