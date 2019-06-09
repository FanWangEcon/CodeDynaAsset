Back to [Fan](https://fanwangecon.github.io)'s
[Dynamic Assets Repository](https://fanwangecon.github.io/CodeDynaAsset/) Table of Content.

# Borrowing and Savings

For dynamic asset models, sometimes we have a dynamic savings choice problem for households and within period static borrowing problem for the firm. We could allow for dynamic borrowing as well as savings.

Generally, the problem becomes a little bit trickier when we have dynamic borrowing as well as savings. One need to consider several small issues carefully.

- *Minimum Income*: Ymin
- *Minimum Consumption*: Cmin
- *Value assigned to invalid choices*: Vinvalid
- *Borrowing bound*: Bbar

## The AZ Problem
> Files for this section are in the [/m_az/](https://github.com/FanWangEcon/CodeDynaAsset/tree/master/m_az) folder in [Fan](https://fanwangecon.github.io/)'s [CodeDynaAsset](https://github.com/FanWangEcon/CodeDynaAsset) repository.

### Common-Grid

**Savings Only**
In the basic AZ problem, Ymin = Zmin*W. Cmin is not needed. Vinvalid is needed to replace invalid choices when we use common-level-grid. Bbar = 0.

The lowest possible utility is reached when households face Ymin and a=0.

- *Ymin*: not parameter, Zmin*W
- *Cmin*: not parameter, 0<Cmin
- *Vinvalid*: parameter, when c < 0 for some a' choices on the common-grid
- *Bbar*: not parameter, constant, a' >= 0

**Save/Borrow + No Default**

If the borrowing problem does not allow for default, the borrowing level is limited by households' ability to repay in future periods unless Bbar is more binding than than.

What at the feasible cash-on-hand levels that could be reached by the choice grid?

The borrowing condition is: borr*(1+r_borr) + Ymin > Bbar. This condition means that a borrowing choice is only valid if when faced with the worst shock tomorrow, the household is still able to repay the debt with new borrowing. This means that: borr_min = (-Ymin/r_borr). This is the *natural borrowing constraint*.

With the Common-Grid choice structure, if we allow the minimum asset point to be = (-Ymin/r_borr), at that point, all choices lead to c = 0. If a = (-Ymin/r_borr), coh = Ymin + (-Ymin/r_borr)(1+r_borr) = -Ymin/r_borr. This means we have to borrow up to the max to get to c = 0. But c = 0, u(c) still not defined.

For all state space points higher than (-Ymin/r_borr), there is a point in the choice grid where c != 0, so with the grid based maximization problem, there is always a choice point where utility is defined. We set all other values to *Vinvalid*, but now every single choice point is undefined. The inclusion of this point is not possible. A solution seems to be to allow for cmin, but that should really only be added if we are allowing for default.

So as a practical matter, once borrowing is allowed, the choice grid should be constructed to allow for points in the savings grid region, zero, and also negative borrowing points, but the points should be just higher than (-Ymin/r_borr) with some threshold utility allowed. perhaps set to the same value as what would be achieved under cmin. but the concept would not be to allow for default, but to approximate the idea of just borrowing a little bit higher than the largest amount of borrowing allowed.

Note that unlike in the risky + safe asset choice problem, here we do not have to worry about issues related to absorbing states. Even when we start at the worst/lowest state level today, the exogenous income process is the same in the next period, so if we draw a good shock, we will exit.

- *Ymin*: not parameter, Zmin*W
- *Cmin*: not parameter, 0<Cmin, although could pick min borrow
- *Vinvalid*: parameter, when c < 0 for some a' choices on the common-grid
    + there is a valid choice point for each state, so Vinvalid does not enter V, it is a algorithm-parameter.
- *Bbar*: parameter, max(ExoBbar, -Ymin/r_borr)


**Save/Borrow + Default**

If default is allowed, what does that mean? That means, even if your current total debt is beyond what you can repay, you will get some minimum level of consumption. Afterwards, you debt is wiped out.

Allowing for default leads to several issues:
1. there must be a ExoBbar (and implicitly a ExoSaveLimit): with cmin, the optimal choice, is to borrow up to negative infinity or save up to infinity. When borrowing up to infinity, capped by consumption floor today, gain from utility tomorrow. When savings up to infinity, capped by consumption floor in the future, gain from utility today. Additionally, need to set Borrowing bounds that are small enough that they are not binding, so they don't actually determine optimal choices.
2. in the standard set-up, all debts are wiped out regardless of how much is borrowed, as long as households are unable to repay. Does this lead to households always borrowing up to the bound? This is not the case, although the borrowing bound if set too high, could become binding as discussed in point 1. The idea is, given that there are different realizations of shocks in the future, the household when borrowing considers in how many states of the world tomorrow it will have to face u(cmin). If the household borrows more, it will face u(cmin) which is not negative infinity, but still not good, in more states. So while it is true that within state, given default the household will not pay a higher cost given cmin is fixed. But ex-ante to the realization of shocks, when borrowing decisions are made, households will be considerate in how much they borrow taking into account in how many states of the world they are willing to take the u(cmin) default utility.

We will clearly still have Ymin, we now need a cmin, which becomes a very important parameter.

How important is cmin in expanding the borrowing choice set?
- the models here are infinite horizon, so not much.
- in a finite horizon model/life-cycle for example, the inclusion of cmin very dramatically changes how much could be borrowed. Previous borrowing bound by total Ymin income from finite future periods. Not that bound no longer exists.

We now set up the choice grid to be between ExoBbar at the lower and and some upper end values. ExoBbar, again should be not binding. If it is binding, it is because it is too high, need to lower it.

One way to set the parameters is that, to avoid the issue of absorbing state, unlike before, where we set borrowing limit by Ymin, we can set it by Ymax. That is, borrowing is find as long as there exists one state of the world tomorrow hen utility is not negative infinity. This is important, Because if at all future states, there is default, that means this will be the case not just tomorrow, but also the day after, and forever. This means that the household has no incentive to limit borrowing. Optimally, the household will borrow up to infinity today. As long as one state in the future remains not in default, there is a tradeoff between today and tomorrow.

This means:
borr*(1+r_borr) + Ymax > BDefaultBar

We have a very similar situation as before, borrowing must be below this point, We don't have to worry about the situation we had before, finding just poitn slighlyt below. This choice grid is a non-binidng point. Households should not be willing to choose this point. THis is a bound because we don't want households to borrow more than this. If they are at this point, their optimal choice is to borrow infinitly. But we are not allowing that, so borrowing up to this poitn is not optimal. Borrowing just a little bit less should be optimal because it will get at least one future states out of u(cmin) which should have very high marginal returns.

Vinvalid becomes tricky now.

- *Ymin*: not parameter, Zmin*W
- *Cmin*: parameter, 0<Cmin
- *Vinvalid*: parameter, when c < 0 for some a' choices on the common-grid
    + for the default states, all choice points are invalid.   
- *Bbar*: parameter, max(BDefaultBar, Bbar)


Distinguishing feature of default vs not-default is under not-default there remains one choice grid point where there is valid utility today, and one could choose a' = 0. Under default, all points invalid. The question is, what is the optimal choice next period and what is the utilitytoday under these default states?
