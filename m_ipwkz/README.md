# The IPWKZ model

- *i*: interpolation
- *p*: percentage choice grids
- *w*: choose aggregate savings w = k' + b'
- *k*: choose k' given w and z
- *z*: some possibly persistent shock process

# Why IPWKZ?

Already have very fast iwkz, what does ipwkz add? Partly it's for finding the proper endogenous asset distributions. With the iwkz model, there is an issue with minimum choice grid level. At the lowest level of cash-on-hand, coh is below the lowest risky investment choice grid point, hence households are stuck there without the ability to generate income. Unlike in the *bz* model, in the benchmark *wkz* model, the shock only has an impact when risky investment is non-zero. This is a problem because at other state space points, households are optimally choosing assets that could land them in this absorbing state. When the model is simulated for distribution of assets, if simulation goes long enough, there is leakage, leading to all mass accumulated at these worst states.

The above issue is happening not because there is an issue with the model, but because there is an issue with the computational approximation of the model. By solving for percentage optimal asset choices given coh and shock, at all levels of coh, households could potentially choose non-zero levels of risky investment and hence rise out of abject poverty. There are no absorbing states. 
