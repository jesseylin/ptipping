### Types of tipping

#### B-tipping

- Bifurcation-induced
- When an external input pushes the system through a bifurcation and the base state either disappears or becomes unstable

#### R-tipping

- Rate-induced
- When external input varies too fast, deviation from the base state crosses a tipping threshold and the dynamics fundamentally change 
- Doesn't require bifurcation

#### N-tipping

- Noise induced
- Random fluctuations push the system past the base state into the domain of attraction of another state
	- 

### The models

Two ecological models for predator-prey dynamics

$N$ is prey and $P$ is predator

### RMA

$\dot{N} = r(t)N(1-\frac{c}{r(t)}N)(\frac{N-\mu}{v+N})-\frac{\alpha NP}{\beta+N}$
$\dot{P} = \chi \frac{\alpha NP}{\beta+N}-\delta P$

- Modified version of Rosenzweig-MacArthur Model
- logistic growth of the prey
- Holling Type-II functional response
	- **Introduction:** In the type II functional response, the rate of prey consumption by a predator rises as prey density increases, but eventually levels off at a plateau (or asymptote) at which the rate of consumption remains constant regardless of increases in prey density (see also  [TYPE I](https://www.nimbios.org/~gross/bioed/bealsmodules/functional1.html)  and TYPE III FUNCTIONAL RESPONSE). Holling's disk equation, named after experiments he performed in which a blindfolded assistant picked up sandpaper disks from a table, describes the type II functional response.
- linear mortality of the predator



Notes from the paper:
- $-\frac{r(t)\mu}{v}$ is low density (negative prey growh rate)
- $\frac{c\mu}{v}$ is "nonlinear modification of low density prey growth)
- $\frac{N-\mu}{v+N}$ accounts for a "strong Allee effect that accounts for negative prey growth rate at low prey population density"
- $\alpha$ predator kill saturation rate
- $\beta$ is half saturation
- $\frac{r(t)}{c}$ is the carrying capacity (max prey population without predators)
- $\chi$ prey to predator conversation rate
- $\delta$ predatory mortality rate

The main differences? Time dependent rate and Allee effect

$r(t)$ is introduced to approximate a varying climate.

Equilibria
- $e_0=(0,0)$, extinction
- $e_1(r) = (\frac{r}{c},0)$, prey only equilibrium
- $e_2 = (\mu,0)$, Allee equilibrium 
- $e_3(r) = (N_3, P_3(r))$
	- where
	- $N_3 = \frac{\delta \beta}{\chi^{\alpha - \delta}}\geq 0$
	- $P_3(r) = \frac{r}{\alpha}(1-\frac{c}{r}N_3)\frac{(\beta+N_3)(N_3-\mu)}{v+N_3}\geq 0$
	
### May

$\dot{N} = r(t)N(1-\frac{c}{r(t)}N)(\frac{N-\mu}{v+N})-\frac{\alpha NP}{\beta+N}$
$\dot{P} = sP(1-\frac{qP}{N+\epsilon})$

- Prey dynamics the same 
- Predatory modified such that 
	- $s$ is low density predator growth rate
	- $\epsilon$ is "other prey"
	- $q$ is minimum prey-to-predator biomass ratio for predatory population growth.
	- $\epsilon/q$ maximum density in absence of $N$

Equilibria
- $e_0=(0,\frac{\epsilon}{q})$, extinction of prey
- $e_1(r) = (\frac{r}{c},0)$, prey only equilibrium
- $e_2 = (\mu,0)$, Allee equilibrium 
- $e_3(r) = (N_3(r), P_3(r))$
- $e_4(r) = (N_3(r), P_3(r))$
	- where
	- $N_3(r)$ and $N_4(r)$ are the nonnegative roots of
		- $N^3-(\mu-\beta+\frac{r}{c}-\frac{\alpha}{cq})N^2-(\beta\mu+\frac{r(\beta-\mu)}{c}-\frac{\alpha(v+\epsilon)}{cq})N+(\frac{r\beta \mu}{c}+\frac{\alpha v \epsilon}{cq})=0$
	- And 
		- $P_3(r)= \frac{N_3(r)+\epsilon}{q}$
		- $P_4(r)= \frac{N_4(r)+\epsilon}{q}$

### Bifurcations

For RMA
- As $r$ increased, the coexistence equilibrium $e_3(r)$ undergoes a supercritical Hopf bifurcation. We get a stable limit cycle $\Gamma(r)$
- At $r_h$,  the lower boundary of the bifurcation $\Gamma_{\min}$ intersects $N = 0$
	- Here, there's a heteroclinic bifurcation $h$, past which the only attractor is $e_0$
- Bifurcation tipping can happen from all phases but it's proportional to where the system spends the most time
- Tipping can't occur when $r(t)$ decreases because there's no bifurcations it can hit

For May
- Also experiences a supercritical Hopf bifurcation along $e_3(r)$ at $H_1$ though as $r$ increases it hits a reverse supercritical Hopf bifurcation at $H_2$ and stops cycling.

[Show the limit cycles for both]

### Phase
- To denote where along the limit cycle we are, we define 
$\phi_\gamma = \tan^{-1}(10^3 \frac{P_\gamma-P_3}{N_\gamma-N_3}) \in [0,2\pi)$
Where $(N_3, P_3)$ is the oscillatory coexistence equilibrium. $10^3$ is included as a factor because $N$ is often on an order of 1000 times larger than $P$ (you need a lot of rabbits to feed a lynx!)
This can be used as a metric anywhere in the neighborhood of $\Gamma(r)$

### $r(t)$

- Added to these models to tune the prey birth rate and carrying capacity.
- Oscillatory functions generally used in the literature, but don't we see a lot of abrupt, or piecewise changes now?
	- Rapid temperature change, loss of habitat, invasive species, etc.
- How the authors generated $r(t)$?
	- Choose magnitude uniformly from $[r_1,r_2]$.
	- Use hypergeometric distribution to determine number of constant years
		- $g(l) = (1-\rho)^l \rho$
		- $\rho = 0.2$ per "actual climate records from four locations in the boreal and deciduous-boreal forest in North America"
		- Defined above average years at Type H and below as Type-L

### Not just B-tipping

- Monte Carlo simulation of RMA system using parameters fit on Lynx-Hare dynamics
- Start at $(N_0,P_0) = (3,0.002)$ to be inside the basin of attraction for $\Gamma(r)$
- Main outputs?
	- Statistics of tipping, ie $r_{\text{pre}}$ and $r_{\text{post}}$
	- Phase/coordinate at which the rate changes occurred
- Primary finding - you'd expect $r_{\text{pre}} <r_h$ and $r_{\text{post}} >r_h$ but this is not always true.
	- Sometimes it is, and the likelihood of B-tipping is proportional to invariant measure
- There's many simulations in the high predator low prey portion of the model that go to extinction but never hit $r=r_h$
- The culprit?

### P-tipping

- Phase-tipping
- Essentially any of the other forms of tipping but selective to certain phases of periodic states 
- There's clear intuition here
	- swinging a rock on a rope and letting go - down or up?
	- An unstable limit cycle, moving from inside to outside quickly

- Here
	- $r(t)$ decreases to cause tipping and does not cross any bifurcations of $\Gamma(r)$
	- When $r_1 <r_h$ it is still observed so really no bifurcations
