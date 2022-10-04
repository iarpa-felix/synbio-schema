DROP TABLE IF EXISTS parts_sequences_plus;
create table parts_sequences_plus as
select *
from parts_sequences;
update parts_sequences_plus
set "sequence" = replace(replace("sequence", ' ', ''), '\n', '');
ALTER TABLE parts_sequences_plus
    ADD seq_len int;
update parts_sequences_plus
set seq_len = length("sequence");
