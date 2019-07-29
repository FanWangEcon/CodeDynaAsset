- 2019-07-27 20:49
  + update parameter testing for models with zr, only has single interest rate for borrowing now. 

## 2019-07-25 08:16

update abz:
* [x] re-run solu, less shock points
* [] re-run value function testers, did not have shocks for interest

- 2019-07-19 15:03
  + cl_mt_pol_coh to cl_mt_coh for state coh unless actually pol_coh based on choices.

# Issues 2019-07-14

- [] 2019-07-14 11:36
  + below zero -20 consumption for ikwpkz + fibs, happens to be cmin/default proportion as well.

# To do 2019-07-11

- [] 2019-07-11 14:24
  + wkz model, when defaulting what are next period asset levels, now w = 0, what about a' and b', Now a' and b' differs depending on informal interest rate because at w=0, optimal choices differs depending on the informal interest rate as well as the productivity shock.

# To do 2019-07-10

- [] 2019-07-10 11:54
  + legends for az and akz graphs, show different legends incorporating r shock if r shock exists.

- [x] 2019-07-10 15:02
  + 2019-07-10 16:03
  + Graph imaginery warning
    > Warning: Using only the real component of complex data.
    > In getRealData (line 52)
    > In scatter (line 56)
    > In ff_az_vf_post_graph (line 257)
    > In ff_az_vf_post (line 155)
    > In ff_abz_vf_vecsv (line 401)
    > In ff_abz_ds_wrapper (line 111)

- [x] 2019-07-10 15:07
  + completed: 2019-07-10 16:03
  + mkdir exists warning, use "if ~exist(dirpath,'dir') mkdir(dirpath); end"

- [ ] 2019-07-10 21:17
  + dimension correction in descriptions of ipwkbz: *(I^k x I^w x M^r) by (M^z)* to *(P^k x I^w x M^r) by (M^z)*
