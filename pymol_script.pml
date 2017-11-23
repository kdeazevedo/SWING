load test_data/2za4.pdb, 2za4
load test_data/1AY7_l_sep.pdb, 1AY7_l_sep
load test_data/1AY7_r.pdb, 1AY7_r
select 2za4, 2za4 and chain A
select 2za4, 2za4 and chain B
cealign 1AY7_r, 2za4
cealign 2za4, 1AY7_l_sep
select tosave, 1AY7_l_sep
save test_data/1AY7_l_sep_aligned.pdb, tosave
quit