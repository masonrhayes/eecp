### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ fc9e304e-850c-11ec-24dd-3bf8c6e392a3
using Distributions, PlutoUI, Symbolics, Latexify, DataFrames, Plots, PlotThemes, ImageShow, ImageIO, FileIO

# ╔═╡ c3a34ad8-3405-4e90-8460-f06234f91610
md"""
# [Environmental Preferences and Technological Choices: Is Market Competition Clean or Dirty?](https://scholar.harvard.edu/aghion/publications/environmental-preferences-and-technological-choices-market-competition-clean-or)

by Aghion, Bénabou, Martin, and Roulet (2021) 

*Presentation by* Mason Hayes and Martí Puig
"""

# ╔═╡ 0f4735b1-f2cd-499c-b995-166872a3a1ef
md"""
#### Follow along by scanning this QR code:
"""

# ╔═╡ 96e82adf-2a20-4949-963e-525a5139aa8e
begin 
	url = download("https://masonrhayes.keybase.pub/projects/pluto/aghion2021_qrcode.png")
	qr = load(url)
end

# ╔═╡ d21e9e92-f812-43d7-890e-ec140bff80ed
PlutoUI.TableOfContents(aside = false)

# ╔═╡ b3938d08-cfcb-4021-a847-3fe6fbd76cbb
md"""
## Abstract

> We investigate the effects of consumers’ environmental concerns and market competition on firms’ decisions to innovate in “clean” technologies. Agents care about their consumption and environmental footprint; firms pursue greener products to soften price competition. Acting as complements, these forces determine R&D, pollution, and welfare. We test the theory using panel data on patents by 8,562 automobile-sector firms in 41 countries, environmental willingness-topay, and competition. As predicted, **exposure to prosocial attitudes fosters clean innovation, all the more so where competition is strong**. Plausible increases in both together can spur it as much as a large fuel-price increase.
"""

# ╔═╡ 851053a7-c31b-4722-a443-a8e212516359
md"""
## Exploring the theoretical model
"""

# ╔═╡ e5f8040c-ae07-47a9-b580-4e1ab0465f66
@variables t y k_f γ δ c z κ # declare some variables

# ╔═╡ 1886e820-9d9f-4a37-9fbc-06c67f055583
md"""
- time t, output y, cumulative investments $k_f$
- consumer environmental preferences $δ$
- size of a leading-edge environmental innovation $γ$
- marginal cost of labor $c$
-  $z$ the investment flow per period = the probability of successful innovation
-  $κ$ some measure of ideas/managerial capacity that makes innovation more costly
"""

# ╔═╡ 037fc3bd-2042-4d5a-943c-289a983909c6
md"""
### Variable definitions

**Emissions**: The production or consumption of one unit of good with environmental quality $q$ generates $x = 1/q$ units of polluting emissions.

**Quality** for a firm $f$ that has the cumulative number of investments $k_f ∈ N$ is given by:
"""

# ╔═╡ 1e8cbf03-ae43-43df-82ee-1e7bcdb1663f
begin
	q_f((γ, k_f)) = γ^k_f
	q_f((γ, k_f))
end

# ╔═╡ 13922ce5-7d66-43dd-a813-0735f31c7173
md"""
Firms choose to invest or not in each period. In any period, firms can copy the other firm's *previous* technologies. Let's assume a duopoly with firms A and B.

Then, in any period $t$, the maximum of $| k_A - k_B |$ must be less than or equal to 1. Since, if it was greater than 1, one of the firms could have copied a better technology than the one it has.
"""

# ╔═╡ b0a2b5d0-d168-4339-8177-4299f2480293
md"""
#### How does quality change with a change in $\gamma$ or $k_f$?
"""

# ╔═╡ cde9049d-224b-4312-bc9c-09eac980098e
md"""
We can differentiate this function wrt. to $\gamma$ or $k_f$ to see how changes in preferences or cumulative investments affect quality
"""

# ╔═╡ 6892aa11-e06c-4a7d-aa6f-f66e8bb61b3d
q_f′(a) = Symbolics.derivative(q_f((γ, k_f)), a; simplify = true);

# ╔═╡ bb213912-cda5-4bf5-9253-69d70f226c65
q_f′(k_f)

# ╔═╡ aa9c1bc2-ae3a-4004-8e63-5300b3189e06
q_f′(γ)

# ╔═╡ 9c7bda65-0996-41d4-a26b-06f290823aae
md"""

So we see that 

``\dfrac{\delta q_f}{\delta k_f}= \gamma^{k_{f}} \log\left( \gamma \right)``, and ``\dfrac{\delta q_f}{\delta \gamma} = \gamma^{(k_{f} - 1)} k_f``, which are both greater than 0.

"""

# ╔═╡ c3bd0b54-adef-4c90-a058-2dbb6b014663
md"""
#### Costs

To produce *1 unit of quality-adjusted output*, a firm must use the following units of labor (with wage normalized to 1):

> Note that we assume **labor** is the only input.
"""

# ╔═╡ 701d66fb-55fc-419f-8c76-91fdf4ed7739
begin 
	L((c, γ, k_f, δ)) = γ^(-k_f * δ) * c;
	L′(a) = Symbolics.derivative(L((c, γ, k_f, δ)) , a);
end;

# ╔═╡ 174ded77-b0e1-4b5a-ae5e-f92887abd80b
L((c, γ, k_f, δ))

# ╔═╡ 5061f3ed-4626-4d21-a087-51e2cfae452d
dLdk_f = L′(k_f) < 0
# labor cost per unit of quality-adjusted output is decreasing in cumulative number of clean innovations

# ╔═╡ 80ba815c-fc9a-49b0-9f43-7b72f7f0496d
md"""
### How do consumers value product quantity & quality?
"""

# ╔═╡ 9e0cc16c-b199-4ff7-b087-1bf00754b9b5
md"""
With $y$ being the quantity produced and $q$ the quality of the product, consumers value a quantity-quality combination of that product at $y q^δ$.

Revenue per sector is normalized to 1.
"""

# ╔═╡ ef7c1d16-d266-47a1-9e7e-21c9a4671e97
u_i((y,a)) = y * q_f((γ, (k_f+a) * δ));

# ╔═╡ a9918e58-0eb0-4fd7-a3dc-584c9a92ccf5
md"""
## The choice to invest
"""

# ╔═╡ 60b97c59-7631-4e5e-b5b6-a97c4cee0653
md"""
For any $z \leq 1$, investing $\kappa z^2/2$ units of labor leads to a probability $z$ of inventing a technology that is $\gamma$ times cleaner; and a $1 - z$ probability of no progress. This means that investing z leads to quality:
"""

# ╔═╡ d74fec10-600d-4c5a-8cf3-3e78d926694b
begin
	Invest((z)) = z * q_f((γ, k_f+1)) + (1-z) * q_f((γ, k_f));
	Invest((z))
end

# ╔═╡ f4fca4ce-1fd4-4bfa-8142-0dbd544e7fa1
Invest′(a) = Symbolics.derivative(Invest((z)), a, simplify = true);

# ╔═╡ 84c4350b-8330-4479-9933-d3751d095aaf
Invest′(γ) # quality is increasing in the size of the leading-edge innovation γ

# ╔═╡ 1bb06086-97a6-4e7f-a103-a4acd0a68803
Invest′(k_f)

# ╔═╡ cc1a7094-eab8-443e-bd44-b83914ef1015
md"""
## Consider an oligopoly with firms A and B. 

And consider an *unleveled sector* where, for example, $k_A = k_B + 1$ which means that firm A has invested $\kappa z^2 /2$ units of labor to produce a technology that is $\gamma$ times cleaner. 

Now, firm A has a product that is more valued by consumers.

"""

# ╔═╡ b49dfe41-ca19-4d79-bccd-cb631c25385d
md"""
Since we know that consumers value firm A's good at $y q^δ = y \gamma^{\delta (k_A)}$, that firm A can capture all demand (by assumption), and that the firm has cost $c$, it can engage in *limit pricing* because it has the following quality advantage over firm B:

"""

# ╔═╡ 1b70417c-e88b-4266-9242-fc2998f33906
(q_f((γ, k_f + 1))/q_f((γ, k_f)))^δ

# ╔═╡ f003f65c-840a-4ae8-994e-5d101eca5c55
md"""

> Note that this quality advantage just equals $γ^δ$.

### Limit pricing

The firm with the lowest price/quality ratio gets all demand and, because of competitive pressure, chooses price so that entry is not profitable. Since it has quality advantage $\gamma^\delta$ compared to its competitor, it can engage in *limit pricing* and choose the monopoly price:
"""

# ╔═╡ c5cebadb-f9ae-4760-b52d-ed7d93173ce8
 pᴹ = γ^δ*c 

# ╔═╡ 84810e98-0a4e-424c-be47-7b77c9cf4ccb
md"""
Demand is:
"""

# ╔═╡ 06f077cc-129c-4577-87d1-89acc387ff9e
yᴹ = 1 / pᴹ

# ╔═╡ 67bc75f5-28df-44f6-97e1-944c50be7b35
md"""
Profits are:
"""

# ╔═╡ ae6f27be-c8cb-4abf-a3ba-6e8d84fa2afe
πᴹ = yᴹ*(pᴹ - c)

# ╔═╡ 5c0da089-c726-498b-9387-0e7de160fd4f
dπᴹd(a) = Symbolics.derivative(πᴹ, a);

# ╔═╡ d1c1168a-fb09-4b12-85af-e713c951cd97
dπᴹd(γ) |> simplify 
# profits are increasing in γ, the size of the leading-edge environmental innovation

# ╔═╡ 58be1c37-928d-4dce-8764-4f661d757351
dπᴹd(δ) |> simplify 
# and profits are increasing in consumer preferences for green technologies

# ╔═╡ 8c0de67b-30c8-4726-be85-d7f76a5d37fe
md"""
## Level of competition
"""

# ╔═╡ b326f0e1-60d9-481b-b860-1c6bb901e2e8
md"""

 $\Delta$ is a measure of the level of competition between the two firms. A higher $\Delta$ implies that each firm's profit when in a duopolistic market is approaching zero (as $\Delta \rightarrow 1$, then $\pi_D \rightarrow 0$. And conversely, as $\Delta$ approaches $\frac{1}{2}$, then the two firms become closer to splitting evenly the monopoly profits, which indicates perfect collusion (lack of competition).

"""

# ╔═╡ fa332e36-0e84-4867-ae8b-3682769e9a3f
π_D(Δ) = (1 - Δ) * πᴹ;

# ╔═╡ 685cbb78-1aed-4b0f-9447-dae359f9a556
p(Δ) = c / (1 - 2*π_D(Δ));

# ╔═╡ 2a76d08a-f383-435c-81f9-849e2143e52a
output(Δ) = 1/p(Δ);

# ╔═╡ bce0b40f-1fde-42e1-ab34-b31e1b80846d
md"""
Duopoly **profits** | Δ =
"""

# ╔═╡ cac3312d-f956-48da-b190-13c494b31d45
md"""
Duopoly **price** | Δ:
"""

# ╔═╡ 7feb16f3-2337-47bf-b97b-2fa935c37466
md"""
Duopoly **output** | Δ:
"""

# ╔═╡ e34d470e-4a14-4aca-8b28-1f23d332eaec
md"""
## Escaping competition through clean innovation
"""

# ╔═╡ 8d13949d-8188-48df-98ae-e9638b7f917f
md"""
When the sector is *leveled* -- when $k_A = k_B$ -- only one of the two firms is given the opportunity to innovate in each period. $\forall z$ s.t. $0 \leq z \leq 1$, the cost of innovation is $\kappa z^2 /2$ and leads to monopoly profits with probability $z$.

So, the firm solves the following maximization problem:
"""

# ╔═╡ 2bdcaa58-fcce-4edb-8cf1-71d388699e7d
md"""
``\max_{z ∈ [0,1]} \{z π^M + (1-z) π_D(Δ) - κz^2 /2\}``
"""

# ╔═╡ 877190a3-2024-4ef0-bda7-7081d6c145fb
md"""
Let's let Julia solve this for us:
"""

# ╔═╡ 20c63144-8102-4204-8b0a-f3e861a50ceb
md"""
Notice that the optimal z, denoted by $\hat z$, is just the difference in a firm's monopoly profit and duopoly profit divided by $\kappa$, but can never exceed 1.
"""

# ╔═╡ 214c6ce1-fd94-460c-8a46-d150e4cc85be
md"""
Simplifying, we can see that $\hat z$ can be written as:
"""

# ╔═╡ eccbf170-1ee3-49fa-8ad3-4f407e86df1d
I(Δ) = Δ*πᴹ/κ # define investment flow
# the investment flow per period is just z(Δ);

# ╔═╡ b2c82550-c3ed-4804-847f-63cc630bf5f8
md"""
## Aggregate flow of investments per period

The flow of investments per period is just the proportion of sectors where innovation will occur:
"""

# ╔═╡ 755a0794-e085-4419-b9fd-68062367b8da
md"""
### Per-period investment is increasing in Δ and δ
#### and the two forces *complement each other*
"""

# ╔═╡ b2c953cf-eb33-45e5-81a5-2474cde70082
md"""
> Recall that $Δ$ is the measure of the level of competition, and $δ$ is a measure of the strength of consumers' social-responsibility concerns
"""

# ╔═╡ ee924789-7648-494c-9bd6-e749d204cf95
@variables Δ;

# ╔═╡ 9090cf3a-5505-4f0d-8876-b39de08bc5f0
π_D(Δ) |> simplify

# ╔═╡ 076e8d22-9bad-41d1-805a-8eff7a3f8789
p(Δ) |> simplify

# ╔═╡ 85db7a2f-6818-4460-be8f-fe5755f9dbe8
output(Δ) |> simplify

# ╔═╡ c188679c-c689-4b84-82a7-94872c4ff4f2
dπdz = Symbolics.derivative(z*πᴹ + (1-z)*π_D(Δ) - κ*z^2 /2, z)

# ╔═╡ daca5c07-366d-44bb-bfd4-37528000768c
ẑ = Symbolics.solve_for(dπdz ~ 0, z)

# ╔═╡ 39ec2265-f90a-4f30-aa3f-20e0814f7afb
ẑ |> simplify

# ╔═╡ 70940ee4-b6cf-4783-acfd-3277b1252daa
I(Δ) |> expand |> simplify

# ╔═╡ 4376cf9e-5239-4174-92e1-89e081da4f37
I′(a) = Symbolics.derivative(I(Δ), a) # make a function to take the derivative of I(Δ);

# ╔═╡ eed08415-9080-4bd9-af9d-6ed3be107c6a
dIdδ = I′(δ) |> Symbolics.simplify_fractions > 0

# ╔═╡ e67aca04-24e0-4395-82c5-321cdc35d55a
dIdΔ = I′(Δ) |> expand |> Symbolics.simplify_fractions > 0

# ╔═╡ 86373d06-c9bc-465b-8e04-dca5ca7b4706
d²IdδdΔ = Symbolics.derivative(I′(δ), Δ)|> simplify > 0 # = d²I/dΔdδ

# ╔═╡ 9c068dd5-385d-44a6-975a-996873bb35bc
md"""
#### Interpretation

The above equations show that **more competition** as well as **stronger preferences for green technologies** each increase the level of investment in such technologies and therefore the level of innovation.
"""

# ╔═╡ 5791f07a-f964-4b4a-a37c-f7fe6aed7560
md"""
> But how do the two influence the levels of *pollution* and *welfare*?
"""

# ╔═╡ f8f14043-37d4-4f6d-a188-1c85c935f5f2
md"""
## Pollution and welfare
"""

# ╔═╡ 6c01741c-8637-416a-bd57-bdcd03c2154c
md"""
Total emissions normalized by total expenditure are given by:

$[1 - I(Δ)] y(Δ) + I(Δ) y^M / γ$

"""

# ╔═╡ aba22dcd-eaca-4b16-9917-f5eaf75839b0
X(Δ) = (1 - I(Δ))*output(Δ) + I(Δ)*yᴹ/γ

# ╔═╡ d4892a0d-016a-4baa-9edd-afde9e26a000
X′(a) = Symbolics.derivative(X(Δ), a); # the derivative of total emission w.r.t variable a (input by user)

# ╔═╡ 87d08e62-2a75-4faf-a418-2670dda64d05
X′(Δ)

# ╔═╡ 3284a4f5-5be5-4c1a-a72a-fa1feb07a699
X′(δ)

# ╔═╡ d7e6229d-6c8e-4a6c-9365-673d506e771f
md"""
The equations here get quite messy, so an immediate interpretation is not very obvious. Let's see if we can play around with some parameters to see how the level of emmissions changes with a change in competition, for example.
"""

# ╔═╡ 9e652832-665b-4d22-9aad-20c07b457c3b
md"""
### Competition and its effects

By increasing competition in *leveled sectors*, those in which firms are neck-and-neck, competition increases pollution.

> As competition increases, ↑ Δ ⟹ ↓ prices, ↑ output ⟹ ↑ pollution since more goods are consumed.
"""

# ╔═╡ 6b7de976-800b-46d0-9875-76e4f8cf2047
md"""
**But**, the fear of lower profits gives firms a high incentive to invest in R&D, which means that the probability of successful innovation $z$ is higher. When these successful innovations occur, emissions are reduced.
"""

# ╔═╡ 9354c637-9f4e-4381-9350-712237372507
md"""
The way that emissions and welfare change depends on the level of $κ$. 
"""

# ╔═╡ 7a259d76-bfd3-4c8e-9958-de900a052d83
κ₁ = πᴹ

# ╔═╡ 3061902f-7c89-4300-8d50-a836beece8db
κ₂ = 1 - γ^(-δ) * (1 + 1/γ)/2

# ╔═╡ 27e24a66-b8da-4721-8fd0-ec274875463d
κ₂ > κ₁

# ╔═╡ bad7baad-2f48-495e-8b0e-6726fcec5371
md"""
## How do emissions change based on changes in:

-  $κ$, the measure of (lack of) ideas/managerial capacity that makes innovation more costly?
- Marginal cost c?
- Level of competition $Δ$?
- Consumer preferences for green technology $δ$?
- The size of the leading-edge environmental innovation $γ$?
"""

# ╔═╡ ce6f05f4-af9f-45c7-811c-f9fc5d3269ad
md"""
> Remember that $γ > 1$ and $δ > 1 \implies c < γ^δ c$
"""

# ╔═╡ 5a1d7422-b6ac-41d0-8051-70fbf8d59236
md"""
#### Visualizing the model

Marginal cost of production c = $(@bind cost Slider(0.50:0.50:10.0, default = 0.50, show_value = true))

Size of leading-edge innovation γ = $(@bind gamma Slider(1.05:0.05:10.0, default = 1.50, show_value = true))

Consumer preferences for greener technology δ = $(@bind delta Slider(1.0:0.05:10.00, default = 1.6, show_value = true))

The cost of innovating κ = $(@bind kappa Slider([κ₂ - κ₁/2, κ₁, (κ₂ - κ₁/2)/0.50, κ₂, (κ₂ + κ₁/2)/1.05, κ₂ + κ₁/2], show_value = true))

"""

# ╔═╡ 4e873006-8c8c-467b-a8a4-4a8a4e1f13d7
md"""
Adjust the limits of the y axis below ↓:
"""

# ╔═╡ 2eca6601-21fa-47ca-a910-30e0a3e7d263
begin
	ylim1 = @bind a Slider(-5:0, default = -1, show_value = true);
	ylim2 = @bind b Slider(0:5, default = 2, show_value = true);
	ylims = [ylim1, ylim2]
end

# ╔═╡ 0d384486-b3ec-4526-8804-d42687644c6d
Delta = @bind Delta Slider(0.5:0.01:0.99, default = 0.75, show_value = true)

# ╔═╡ e0a1f62a-de89-4247-b15d-c7c88a804988
md"""
By tweaking the model, we can see the main results of the paper:

##### Proposition 2

- for $κ < κ_2− κ_1/2$, aggregate pollution $X(∆)$ decreases monotonically;
- for $κ > κ_2− κ_1/2$, $X(∆)$ increases monotonically
- for $κ ∈ (κ_2 − κ_1/2, κ_2 +κ_1/2)$, $X(∆)$ is hump-shaped; moreover, it is minimized at $∆ = 1$ (versus $∆ = 1/2$) if and only if $κ < κ_2$;
- for all $κ ∈ [κ_1, κ_2]$, $X(∆)$ is minimized at $∆ = 1$.

***
##### Proposition 3
- Aggregate pollution $X(∆)$ decreases with consumer’s social-responsibility concern $δ$.

***

##### Proposition 4
- For $κ ∈ [κ_1, κ_2 −κ_1/2]$, social welfare $W$ increases monotonically with competition $Δ$.
"""

# ╔═╡ e0bd7688-f1ed-44e5-ab1b-42631c2a21b7
md"""
## Social welfare
"""

# ╔═╡ c18a134c-e164-4a48-9ec5-553d17b53eb2
U = (1 - I(Δ))*log(output(Δ)) + I(Δ)*log(γ^δ * yᴹ)

# ╔═╡ 4664d0b6-cff9-4bf3-801f-9ade3e9de2ac
begin 
	emissions = round((substitute(substitute(X(Delta), Dict(κ => kappa)), Dict((γ => gamma, δ => delta, c => cost))) |> simplify).val, digits = 4);
	welfare = round((substitute(substitute(U, Dict((κ => kappa, Δ => Delta))), Dict((γ => gamma, δ => delta, c => cost))) |> simplify).val, digits = 4);
end;

# ╔═╡ 40c38155-60b5-416e-a39c-39c30b974222
md"""
Assuming that ``δ =`` $delta, ``γ =`` $gamma, ``c =`` $cost, ``Δ =`` $Delta, and ``κ =`` $kappa, then:

#### Emissions $(emissions)

#### Welfare: $(welfare)

"""

# ╔═╡ f78d9a27-5cc5-4846-95f0-057252eaa762
begin 
	emissions_v_competition(x) = round((substitute(substitute(X(x), Dict(κ => kappa)), Dict((γ => gamma, δ => delta, c => cost))) |> simplify).val, digits = 4);
	welfare_v_competition(x) = round((substitute(substitute(U, Dict((κ => kappa, Δ => x))), Dict((γ => gamma, δ => delta, c => cost))) |> simplify).val, digits = 4);
end;

# ╔═╡ 65c592fe-cf9a-4c18-b779-b719af9e5220
my_data = DataFrame(emissions = [emissions_v_competition(i) for i in 0.50:0.01:1.0],
	welfare = [welfare_v_competition(i) for i in 0.50:0.01:1.0],
	Δ = 0.50:0.01:1.0);

# ╔═╡ a295ec9e-6a00-4b9e-8254-9ae70cf2cf51
begin
	plot(my_data.Δ, my_data.emissions, label = "Emissions", color = "white", theme(:juno))
	plot!(my_data.Δ, my_data.welfare, label = "Welfare", color = "green")
	plot!(legend = :topright)
	title!("Pollution & Welfare vs Competition")
	xlabel!("Competition (Δ)")
	ylims!(a, b)
	annotate!([(0.9584, a, ("More competition →", 8, :bottom, :juno))])
	annotate!([(0.546, a, ("← Less competition", 8, :bottom, :juno))])
	annotate!([(0.57, b - 0.10, ("δ = $delta, γ = $gamma, c = $cost", 8, :top, :juno))])
	
end

# ╔═╡ ca90d22c-b027-428c-9571-16bcd65ce05a
Symbolics.derivative(U, Δ)

# ╔═╡ 25c5a5d7-6ede-472f-81db-fa9daf1600f1
md"""
## Empirical Analysis

Does the model hold in reality?

In short -- yes! Let's look at a quick overview of the empirical findings from the paper:

"""

# ╔═╡ b2fc04ae-76e9-4435-bb0e-543141dbbb56
html"""
      <iframe width="560" height="315" src="https://siasky.net/AAChA3GA7Rk9jo8MqiHmJUohJE1omIkwNVAjX9i5-zrz1A" frameborder="0" allowfullscreen></iframe>
      </iframe>
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
ImageIO = "82e4d734-157c-48bb-816b-45c225c6df19"
ImageShow = "4e3cecfd-b093-5904-9786-8bbb286a6a31"
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
PlotThemes = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"

[compat]
DataFrames = "~1.3.2"
Distributions = "~0.25.46"
FileIO = "~1.13.0"
ImageIO = "~0.6.1"
ImageShow = "~0.3.3"
Latexify = "~0.15.9"
PlotThemes = "~2.0.1"
Plots = "~1.25.9"
PlutoUI = "~0.7.34"
Symbolics = "~4.3.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "6f1d9bc1c08f9f4a8fa92e3ea3cb50153a1b40d4"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.1.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.AbstractTrees]]
git-tree-sha1 = "03e0550477d86222521d254b741d470ba17ea0b5"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.3.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgCheck]]
git-tree-sha1 = "a3a402a35a2f7e0b87828ccabbd5ebfbebe356b4"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.3.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "1ee88c4c76caa995a885dc2f22a5d548dfbbc0ba"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "3.2.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AutoHashEquals]]
git-tree-sha1 = "45bb6705d93be619b81451bb2006b7ee5d4e4453"
uuid = "15f4f7f2-30c1-5605-9d31-71845cf9641f"
version = "0.2.0"

[[deps.BangBang]]
deps = ["Compat", "ConstructionBase", "Future", "InitialValues", "LinearAlgebra", "Requires", "Setfield", "Tables", "ZygoteRules"]
git-tree-sha1 = "d648adb5e01b77358511fb95ea2e4d384109fac9"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.3.35"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Baselet]]
git-tree-sha1 = "aebf55e6d7795e02ca500a689d326ac979aaf89e"
uuid = "9718e550-a3fa-408a-8086-8db961cd8217"
version = "0.1.1"

[[deps.Bijections]]
git-tree-sha1 = "705e7822597b432ebe152baa844b49f8026df090"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.3"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "215a9aa4a1f23fbd05b92769fdd62559488d70e9"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.1"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "f9982ef575e19b0e5c7a98c6e75ee496c0f73a93"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.12.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "12fc73e5e0af68ad3137b886e3f7c1eacfca2640"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.17.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "3f1f500312161f1ae067abe07d13b40f78f32e07"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.8"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "68a0743f578349ada8bc911a5cbd5a2ef6ed6d1f"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.0"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "44c37b4636bc54afac5c574d2d02b625349d6582"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.41.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.CompositeTypes]]
git-tree-sha1 = "d5b014b216dc891e81fea299638e4c10c657b582"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.2"

[[deps.CompositionsBase]]
git-tree-sha1 = "455419f7e328a1a2493cabc6428d79e951349769"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.1"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "ae02104e835f219b8930c7664b8012c93475c340"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.3.2"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DefineSingletons]]
git-tree-sha1 = "0fba8b706d0178b4dc7fd44a96a92382c9065c2c"
uuid = "244e2a9f-e319-4986-a169-4d1fe445cd52"
version = "0.1.2"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "84083a5136b6abf426174a58325ffd159dd6d94f"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.9.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "2e97190dfd4382499a4ac349e8d316491c9db341"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.46"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "StaticArrays", "Statistics"]
git-tree-sha1 = "5f5f0b750ac576bcf2ab1d7782959894b304923e"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.5.9"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.DynamicPolynomials]]
deps = ["DataStructures", "Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "74e63cbb0fda19eb0e69fbe622447f1100cd8690"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.4.3"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.EllipsisNotation]]
deps = ["ArrayInterface"]
git-tree-sha1 = "d7ab55febfd0907b285fbf8dc0c73c0825d9d6aa"
uuid = "da5c29d0-fa7d-589e-88eb-ea29b0a81949"
version = "1.3.0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ae13fcbc7ab8f16b0856729b050ef0c446aa3492"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.4+0"

[[deps.ExprTools]]
git-tree-sha1 = "56559bbef6ca5ea0c0818fa5c90320398a6fbf8d"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.8"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "80ced645013a5dbdc52cf70329399c35ce007fae"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.13.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "deed294cde3de20ae0b2e0355a6c4e1c6a5ceffc"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.12.8"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "51d2dfe8e590fbd74e7a842cf6d13d8a2f45dc01"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.6+0"

[[deps.GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "4a740db447aae0fbeb3ee730de1afbb14ac798a1"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.63.1"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "aa22e1ee9e722f1da183eb33370df4c1aeb6c2cd"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.63.1+0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "1c5a84319923bea76fa145d49e93aa4394c73fc2"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.1"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "b51bb8cae22c66d0f6357e3bcb6363145ef20835"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.5"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "9a5c62f231e5bba35695a20988fc7cd6de7eeb5a"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.3"

[[deps.ImageIO]]
deps = ["FileIO", "JpegTurbo", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "464bdef044df52e6436f8c018bea2d48c40bb27b"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.1"

[[deps.ImageShow]]
deps = ["Base64", "FileIO", "ImageBase", "ImageCore", "OffsetArrays", "StackViews"]
git-tree-sha1 = "d0ac64c9bee0aed6fdbb2bc0e5dfa9a3a78e3acc"
uuid = "4e3cecfd-b093-5904-9786-8bbb286a6a31"
version = "0.3.3"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "87f7662e03a649cffa2e05bf19c303e168732d3e"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.2+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[deps.IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[deps.InitialValues]]
git-tree-sha1 = "4da0f88e9a39111c2fa3add390ab15f3a44f3ca3"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.3.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IntervalSets]]
deps = ["Dates", "EllipsisNotation", "Statistics"]
git-tree-sha1 = "3cc368af3f110a767ac786560045dceddfc16758"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.5.3"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[deps.InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "a77b273f1ddec645d1b7c4fd5fb98c8f90ad10a5"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.1"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LabelledArrays]]
deps = ["ArrayInterface", "ChainRulesCore", "LinearAlgebra", "MacroTools", "StaticArrays"]
git-tree-sha1 = "3696fdc1d3ef6e4d19551c92626066702a5db91c"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.7.1"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a8f4f279b6fa3c3c4f1adadd78a621b13a506bce"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.9"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "340e257aada13f95f98ee352d316c3bed37c8ab9"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "e5718a00af0ab9756305a0392832c8952c7426c1"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.6"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.MappedArrays]]
git-tree-sha1 = "e8b359ef06ec72e8c030463fe02efe5527ee5142"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.1"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[deps.Metatheory]]
deps = ["AutoHashEquals", "DataStructures", "Dates", "DocStringExtensions", "Parameters", "Reexport", "TermInterface", "ThreadsX", "TimerOutputs"]
git-tree-sha1 = "0886d229caaa09e9f56bcf1991470bd49758a69f"
uuid = "e9d8d322-4543-424a-9be4-0cc815abe26c"
version = "1.3.3"

[[deps.MicroCollections]]
deps = ["BangBang", "InitialValues", "Setfield"]
git-tree-sha1 = "6bb7786e4f24d44b4e29df03c69add1b63d88f01"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.1.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "b34e3bc3ca7c94914418637cb10cc4d1d80d877d"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.3"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.MultivariatePolynomials]]
deps = ["DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "fa6ce8c91445e7cd54de662064090b14b1089a6d"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.4.2"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "73deac2cbae0820f43971fad6c08f6c4f2784ff2"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "0.3.2"

[[deps.NaNMath]]
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore"]
git-tree-sha1 = "18efc06f6ec36a8b801b23f076e3c6ac7c3bf153"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "043017e0bdeff61cfbb7afeb558ab29536bbb5ed"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.8"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "923319661e9a22712f24596ce81c54fc0366f304"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.1.1+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "648107615c15d4e09f7eca16307bc821c1f718d8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.13+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "ee26b350276c51697c9c2d88a072b339f9f03d73"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.5"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "eb4dbb8139f6125471aa3da98fb70f02dc58e49c"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.3.14"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "03a7a85b76381a3d04c7a1656039197e70eda03d"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.11"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "0b5cfbb704034b5b4c1869e36634438a047df065"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.1"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "a7a7e1a88853564e551e4eba8650f8c38df79b37"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.1.1"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "6f1b25e8ea06279b5689263cc538f51331d7ca17"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.1.3"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "1d0a11654dbde41dc437d6733b68ce4b28fbe866"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.25.9"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "8979e9802b4ac3d58c503a20f2824ad67f9074dd"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.34"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "db3a23166af8aebf4db5ef87ac5b00d36eb771e2"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "2cf929d64681236a2e074ffafb8d568733d2e6af"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.3"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "afadeba63d90ff223a6a48d2009434ecee2ec9e8"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.1"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "37c1631cb3cc36a535105e6d5557864c82cd8c2b"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.5.0"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "ChainRulesCore", "DocStringExtensions", "FillArrays", "LinearAlgebra", "RecipesBase", "Requires", "StaticArrays", "Statistics", "ZygoteRules"]
git-tree-sha1 = "5144e1eafb2ecc75765888a4bdcd3a30a6a08b14"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.24.1"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Referenceables]]
deps = ["Adapt"]
git-tree-sha1 = "e681d3bfa49cd46c3c161505caddf20f0e62aaa9"
uuid = "42d2dcc6-99eb-4e98-b66c-637b7d73030e"
version = "0.1.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "cdbd3b1338c72ce29d9584fdbe9e9b70eeb5adca"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.1.3"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "cdc1e4278e91a6ad530770ebb327f9ed83cf10c4"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.3"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.SciMLBase]]
deps = ["ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "RecipesBase", "RecursiveArrayTools", "StaticArrays", "Statistics", "Tables", "TreeViews"]
git-tree-sha1 = "f4862c0cb4e34ed182718221028ba1bf50742108"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.26.1"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "0afd9e6c623e379f593da01f20590bacc26d1d14"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.8.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "8fb59825be681d451c246a795117f317ecbcaa28"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.2"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "8d0c8e3d0ff211d9ff4a0c2307d876c99d10bdf1"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.2"

[[deps.SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "39c9f91521de844bad65049efd4f9223e7ed43f9"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.14"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "7f5a513baec6f122401abfc8e9c074fdac54f6c1"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.4.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "a635a9333989a094bddc9f940c04c549cd66afcf"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.3.4"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
git-tree-sha1 = "d88665adc9bcf45903013af0982e2fd05ae3d0a6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "51383f2d367eb3b444c961d485c565e4c0cf4ba0"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.14"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "f35e1879a71cca95f4826a14cdbf0b9e253ed918"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.15"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "d21f2c564b21a202f4677c0fba5b5ee431058544"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.4"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "Metatheory", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "TermInterface", "TimerOutputs"]
git-tree-sha1 = "bfa211c9543f8c062143f2a48e5bcbb226fd790b"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "0.19.7"

[[deps.Symbolics]]
deps = ["ArrayInterface", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "IfElse", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "Metatheory", "NaNMath", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "TermInterface", "TreeViews"]
git-tree-sha1 = "074e08aea1c745664da5c4b266f50b840e528b1c"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "4.3.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "bb1064c9a84c52e277f1096cf41434b675cd368b"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.TermInterface]]
git-tree-sha1 = "7aa601f12708243987b88d1b453541a75e3d8c7a"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "0.2.3"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadsX]]
deps = ["ArgCheck", "BangBang", "ConstructionBase", "InitialValues", "MicroCollections", "Referenceables", "Setfield", "SplittablesBase", "Transducers"]
git-tree-sha1 = "6dad289fe5fc1d8e907fa855135f85fb03c8fa7a"
uuid = "ac1d9e8a-700a-412c-b207-f0111f4b6c0d"
version = "0.1.9"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "991d34bbff0d9125d93ba15887d6594e8e84b305"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.5.3"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "97e999be94a7147d0609d0b9fc9feca4bf24d76b"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.15"

[[deps.Transducers]]
deps = ["Adapt", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "Requires", "Setfield", "SplittablesBase", "Tables"]
git-tree-sha1 = "1cda71cc967e3ef78aa2593319f6c7379376f752"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.72"

[[deps.TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "66d72dc6fcc86352f01676e8f0f698562e60510f"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.23.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "78736dab31ae7a53540a6b752efc61f77b304c5b"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.8.6+1"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╠═fc9e304e-850c-11ec-24dd-3bf8c6e392a3
# ╟─c3a34ad8-3405-4e90-8460-f06234f91610
# ╟─0f4735b1-f2cd-499c-b995-166872a3a1ef
# ╟─96e82adf-2a20-4949-963e-525a5139aa8e
# ╟─d21e9e92-f812-43d7-890e-ec140bff80ed
# ╟─b3938d08-cfcb-4021-a847-3fe6fbd76cbb
# ╟─851053a7-c31b-4722-a443-a8e212516359
# ╠═e5f8040c-ae07-47a9-b580-4e1ab0465f66
# ╟─1886e820-9d9f-4a37-9fbc-06c67f055583
# ╟─037fc3bd-2042-4d5a-943c-289a983909c6
# ╟─1e8cbf03-ae43-43df-82ee-1e7bcdb1663f
# ╟─13922ce5-7d66-43dd-a813-0735f31c7173
# ╟─b0a2b5d0-d168-4339-8177-4299f2480293
# ╟─cde9049d-224b-4312-bc9c-09eac980098e
# ╠═6892aa11-e06c-4a7d-aa6f-f66e8bb61b3d
# ╠═bb213912-cda5-4bf5-9253-69d70f226c65
# ╠═aa9c1bc2-ae3a-4004-8e63-5300b3189e06
# ╟─9c7bda65-0996-41d4-a26b-06f290823aae
# ╟─c3bd0b54-adef-4c90-a058-2dbb6b014663
# ╠═174ded77-b0e1-4b5a-ae5e-f92887abd80b
# ╟─701d66fb-55fc-419f-8c76-91fdf4ed7739
# ╟─5061f3ed-4626-4d21-a087-51e2cfae452d
# ╟─80ba815c-fc9a-49b0-9f43-7b72f7f0496d
# ╟─9e0cc16c-b199-4ff7-b087-1bf00754b9b5
# ╠═ef7c1d16-d266-47a1-9e7e-21c9a4671e97
# ╟─a9918e58-0eb0-4fd7-a3dc-584c9a92ccf5
# ╟─60b97c59-7631-4e5e-b5b6-a97c4cee0653
# ╟─d74fec10-600d-4c5a-8cf3-3e78d926694b
# ╟─f4fca4ce-1fd4-4bfa-8142-0dbd544e7fa1
# ╠═84c4350b-8330-4479-9933-d3751d095aaf
# ╠═1bb06086-97a6-4e7f-a103-a4acd0a68803
# ╟─cc1a7094-eab8-443e-bd44-b83914ef1015
# ╟─b49dfe41-ca19-4d79-bccd-cb631c25385d
# ╠═1b70417c-e88b-4266-9242-fc2998f33906
# ╟─f003f65c-840a-4ae8-994e-5d101eca5c55
# ╟─c5cebadb-f9ae-4760-b52d-ed7d93173ce8
# ╟─84810e98-0a4e-424c-be47-7b77c9cf4ccb
# ╟─06f077cc-129c-4577-87d1-89acc387ff9e
# ╟─67bc75f5-28df-44f6-97e1-944c50be7b35
# ╠═ae6f27be-c8cb-4abf-a3ba-6e8d84fa2afe
# ╟─5c0da089-c726-498b-9387-0e7de160fd4f
# ╠═d1c1168a-fb09-4b12-85af-e713c951cd97
# ╠═58be1c37-928d-4dce-8764-4f661d757351
# ╟─8c0de67b-30c8-4726-be85-d7f76a5d37fe
# ╟─b326f0e1-60d9-481b-b860-1c6bb901e2e8
# ╠═fa332e36-0e84-4867-ae8b-3682769e9a3f
# ╠═685cbb78-1aed-4b0f-9447-dae359f9a556
# ╠═2a76d08a-f383-435c-81f9-849e2143e52a
# ╟─bce0b40f-1fde-42e1-ab34-b31e1b80846d
# ╟─9090cf3a-5505-4f0d-8876-b39de08bc5f0
# ╟─cac3312d-f956-48da-b190-13c494b31d45
# ╟─076e8d22-9bad-41d1-805a-8eff7a3f8789
# ╟─7feb16f3-2337-47bf-b97b-2fa935c37466
# ╟─85db7a2f-6818-4460-be8f-fe5755f9dbe8
# ╟─e34d470e-4a14-4aca-8b28-1f23d332eaec
# ╟─8d13949d-8188-48df-98ae-e9638b7f917f
# ╟─2bdcaa58-fcce-4edb-8cf1-71d388699e7d
# ╟─877190a3-2024-4ef0-bda7-7081d6c145fb
# ╟─c188679c-c689-4b84-82a7-94872c4ff4f2
# ╟─20c63144-8102-4204-8b0a-f3e861a50ceb
# ╠═daca5c07-366d-44bb-bfd4-37528000768c
# ╟─214c6ce1-fd94-460c-8a46-d150e4cc85be
# ╠═39ec2265-f90a-4f30-aa3f-20e0814f7afb
# ╠═eccbf170-1ee3-49fa-8ad3-4f407e86df1d
# ╟─b2c82550-c3ed-4804-847f-63cc630bf5f8
# ╠═70940ee4-b6cf-4783-acfd-3277b1252daa
# ╠═4376cf9e-5239-4174-92e1-89e081da4f37
# ╟─755a0794-e085-4419-b9fd-68062367b8da
# ╟─b2c953cf-eb33-45e5-81a5-2474cde70082
# ╟─e67aca04-24e0-4395-82c5-321cdc35d55a
# ╟─eed08415-9080-4bd9-af9d-6ed3be107c6a
# ╟─86373d06-c9bc-465b-8e04-dca5ca7b4706
# ╠═ee924789-7648-494c-9bd6-e749d204cf95
# ╟─9c068dd5-385d-44a6-975a-996873bb35bc
# ╟─5791f07a-f964-4b4a-a37c-f7fe6aed7560
# ╟─f8f14043-37d4-4f6d-a188-1c85c935f5f2
# ╟─6c01741c-8637-416a-bd57-bdcd03c2154c
# ╠═aba22dcd-eaca-4b16-9917-f5eaf75839b0
# ╠═d4892a0d-016a-4baa-9edd-afde9e26a000
# ╠═87d08e62-2a75-4faf-a418-2670dda64d05
# ╠═3284a4f5-5be5-4c1a-a72a-fa1feb07a699
# ╟─d7e6229d-6c8e-4a6c-9365-673d506e771f
# ╟─9e652832-665b-4d22-9aad-20c07b457c3b
# ╟─6b7de976-800b-46d0-9875-76e4f8cf2047
# ╟─9354c637-9f4e-4381-9350-712237372507
# ╟─7a259d76-bfd3-4c8e-9958-de900a052d83
# ╟─3061902f-7c89-4300-8d50-a836beece8db
# ╠═27e24a66-b8da-4721-8fd0-ec274875463d
# ╟─bad7baad-2f48-495e-8b0e-6726fcec5371
# ╟─ce6f05f4-af9f-45c7-811c-f9fc5d3269ad
# ╟─5a1d7422-b6ac-41d0-8051-70fbf8d59236
# ╟─a295ec9e-6a00-4b9e-8254-9ae70cf2cf51
# ╟─4e873006-8c8c-467b-a8a4-4a8a4e1f13d7
# ╟─2eca6601-21fa-47ca-a910-30e0a3e7d263
# ╟─0d384486-b3ec-4526-8804-d42687644c6d
# ╟─40c38155-60b5-416e-a39c-39c30b974222
# ╟─4664d0b6-cff9-4bf3-801f-9ade3e9de2ac
# ╟─f78d9a27-5cc5-4846-95f0-057252eaa762
# ╟─65c592fe-cf9a-4c18-b779-b719af9e5220
# ╟─e0a1f62a-de89-4247-b15d-c7c88a804988
# ╟─e0bd7688-f1ed-44e5-ab1b-42631c2a21b7
# ╠═c18a134c-e164-4a48-9ec5-553d17b53eb2
# ╠═ca90d22c-b027-428c-9571-16bcd65ce05a
# ╟─25c5a5d7-6ede-472f-81db-fa9daf1600f1
# ╟─b2fc04ae-76e9-4435-bb0e-543141dbbb56
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002