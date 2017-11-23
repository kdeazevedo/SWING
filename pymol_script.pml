load test_data/2za4.pdb, test_data/2za4
load test_data/1AY7_l_sep.pdb, test_data/1AY7_l_sep
load test_data/1AY7_r.pdb, test_data/1AY7_r
select test_data/2za4_rec, test_data/2za4 and chain A
select test_data/2za4_lig, test_data/2za4 and chain B
cealign test_data/1AY7_r, test_data/2za4_rec
cealign test_data/2za4_lig, test_data/1AY7_l_sep
select tosave, test_data/1AY7_l_sep
save test_data/1AY7_l_sep_aligned.pdb, tosave
quit