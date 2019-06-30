Support function specific to formal and informal choice problems in *A Choice Amongst Many: Household Borrowing in a Setting with Multiple Providers* ([**Robert M. Townsend**](http://www.robertmtownsend.net/) and [**Fan Wang**](https://fanwangecon.github.io/) 2019).

See [Fan](https://fanwangecon.github.io)'s [CodeDynaAsset](https://github.com/FanWangEcon/CodeDynaAsset)'s webpage for [table of contents](https://fanwangecon.github.io/CodeDynaAsset/).

# What do the functions do

- ffs_fibs_min_c_cost.m: determines which combinations of formal and informal options are optimal conditional on bridge choices. 
- ffs_fibs_inf_bridge.m: deals with the bridge problem itself

- ffs_fibs_min_c_cost_bridge.m: overall, which formal and informal are optimal, invokes ffs_fibs_min_c_cost.m, deals with bridge loans

- ffs_for_br_block_match.m: for formal borrowing, given borrowing choice, what is the closest block.
- ffs_for_br_block_gen.m: generate formal borrowing blocks
