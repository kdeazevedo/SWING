load 2za4.pdb, 2za4
load 1AY7_l_sep.pdb, 1AY7_l_sep
load 1AY7_r.pdb, 1AY7_r
select 2za4_rec, 2za4 and chain A
select 2za4_lig, 2za4 and chain B
cealign 1AY7_r, 2za4_rec
cealign 2za4_lig, 1AY7_l_sep
select tosave, 1AY7_l_sep
save 1AY7_l_sep_aligned.pdb, tosave
quit