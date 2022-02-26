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

# ╔═╡ 2edf0470-88f7-11ec-35a5-c77833b3aad4
using DataFrames, StatsBase, PlutoUI, Dates

# ╔═╡ c160c0bc-0b77-4206-a97c-a0da3d662939
html"""
<body>
<section id="energy-economics-climate-policy-homework" class="level2">
<h2>Energy Economics &amp; Climate Policy Homework</h2>
<p>by Mason Ross Hayes on 2022-02-26</p>
<p>done in  <a href="https://julialang.org"> <img src="https://user-images.githubusercontent.com/39578155/142780578-be2c8aa4-359c-43aa-9d3a-655d4b938f7d.png" width="60" height="38" /></a> using <a href="https://plutojl.org">Pluto.jl</a>🎈</p>
</section>
</body>
</html>

"""

# ╔═╡ 5f6503f1-a50c-454a-baca-55a898146f57
md"""

## Solving the exercise

Let's create a dataframe of the data we are given:
"""

# ╔═╡ 88aea7bd-3be8-4a21-8291-8e8d3a3aeca9
df = sort(
	DataFrame(
	θ = ["H", "M", "L"],
	a = [85000.0, 60000.0, 40000.0],
	hours = [300.0, 5700.0, 2760.0],
	fuel_cost = repeat([72.0], 3),
	cap_cost = repeat([6.85], 3)
	), :a)

# ╔═╡ a28c77ef-e445-4ac2-b183-c4da98f3a4c2
md"""
Now let's set the variable price and the parameter b to be flexible:
"""

# ╔═╡ ebda0a22-0d15-4eca-b7c8-72a80e27819b
p = @bind p Slider(0.0:10.0:800.0, default = 100.0, show_value = true)

# ╔═╡ d3818f0c-2223-4e04-b5f4-50e632087f61
md"""
Using the slider below, we can adjust b until the weighted elasticity equals 0.01 when price = 100 euros per MWh:
"""

# ╔═╡ 0deaef71-dc51-4aa1-943d-933070f8f89b
b = @bind b Slider(2.00:0.005:7.00, show_value = true, default = 5.17)

# ╔═╡ 7b23c101-6b8f-4ecf-878b-0c6b44c3302d
## Notice you can use the arrow keys to be more precise
weighted_elasticity = round(sum(b.*p ./ df.D .* df.time_share), digits = 5)

# ╔═╡ ff0d78ba-7597-476c-8549-bfd1577e092c
md"""
## Case θ = $(@bind θ Select([1,2,3], default = 3))

"""

# ╔═╡ d26c71f8-1a0b-49f8-baed-8dd808dfc08e
md"""
There are 3 potential cases: ``θ ∈ [L, M, H]``, corresponding to cases 1, 2, and 3, respectively. For each case, we have the demand constraint that the optimal capacity must be greater than the demand at the peak period.

Let's refer to any ``θ ∈ [L, M, H]`` as a *period*. Then, for any period, the demand in that period must be less than or equal to the total capacity. If this is satisfied in period 3, then it is always satisfied, since period 3 corresponds to ``θ = H``, the period of highest demand.
"""

# ╔═╡ 64bcbe4e-7c16-4959-ad60-d7c64142b4f4
md"""
### The optimal K̂
"""

# ╔═╡ d52ae61d-6279-4da2-97dd-934d0f97d4b2
# the values were obtained by selecting Case θ and then adjusting K in the cell below to solve the equation

K̂ = [54147., 60825., 83593.]

# ╔═╡ b1c42a12-f33c-4c5f-a94d-2867c83b197a
@bind K NumberField(40000.:100.:100000., default = 83593)

# ╔═╡ a14f7ba2-04be-44ad-853c-9c4732d45109
md"""
#### Adjust ``K`` above so that the expected dual price of the capacity equals its fixed cost per hour (6.85€)
"""

# ╔═╡ a53208d1-7a33-4950-8742-d4430fe621ae
md"""
### Functions and further calculations
"""

# ╔═╡ a0186edd-8afa-46db-8b34-19e275214b5c
md"""
Define functions for demand, elasticity, surplus, and price.
"""

# ╔═╡ a7a3a395-a91a-40d2-80da-3e639d3c39b9
begin
	D((a, b, p)) = a .- b.*p # Demand function
	ε_w((b,p,D,s)) = b.*p./D .* s # = weighted elasticity function
	S((a, b, p)) = 0.5 .* (a.^2 .- (b .* p).^2) # gross surplus
	P((a, b, K)) = (a .- K)/b # price function
end;

# ╔═╡ ce53326a-2a8a-4bb2-9b5f-c1ec1c2114a7
md"""
## The solution

- Weights ``f(θ) =`` ($(round(df.time_share[1], digits = 3)), $(round(df.time_share[2], digits = 3)), $(round(df.time_share[3], digits = 3))) for ``θ = {L, M, H}``, respectively.
- ``b^* =`` $b gives weighted elasticity of 0.01 at price of 100€/MWh. 
- The optimal ``K^* \approx`` $(K̂[3]) MW, which corresponds to a price of $(round(P((df.a[3], b, K̂[3])), digits = 3)) €/MWh.
  - With CCGT of power 450 MW, this corresponds to $(ceil(K̂[3]/450)) power stations.
"""

# ╔═╡ 94da671b-ca21-487e-8b80-d43e82bdbb3e
md"""

### When ``K =`` $(K̂[θ]) and ``θ =`` $(θ), then the optimal price is: $(round(P((df.a[θ], b, K̂[θ])), digits = 3))
"""

# ╔═╡ 83c7f56e-defe-4f72-9e51-44bdba2b1a0c
md"""

Is optimal capacity $(K̂[θ]) less than demand $(D((df.a, 5.17, p))[θ])? 

**$(K̂[θ] <= D((df.a, 5.17, p))[θ])**
"""

# ╔═╡ 6c2648a8-e894-4f35-be28-d3f8aab665a1
# ↑ adjust the value of K above ↑ to make this cell equal to 6.85
round(df.time_share[θ:3] .* (P((df.a[θ:3], b, K)) - df.fuel_cost[θ:3]) |> sum, digits = 3)

# ╔═╡ a88188b7-6b29-47cf-a34c-71473509cc79
md"""

## Case 1

When the capacity constraint is binding even under low demand ($θ = L$), then the optimal ``K^* \approx`` $(K̂[1]) MW. But this is higher than $(D((df.a, 5.17, p))[1]), so cannot be optimal.

## Case 2

So, we know we need more capacity than in Case 1. Let's consider case 2, under periods of medium demand (``θ = M``). Then, the optimal ``K^* \approx`` $(K̂[2]).

But again, this is higher than the demand $(D((df.a, 5.17, p))[2]).

## Case 3

So, we know we need more capacity than in Case 2 as well. Now in case 3, under periods of peak demand (``θ = H``). Then, the optimal ``K^* \approx`` $(K̂[3]) which lies within the range of medium demand and peak demand, ``K ∈`` $(D((df.a, 5.17, p))[2], D((df.a, 5.17, p))[3]).

"""

# ╔═╡ 2dc1da6e-3da5-4681-b7bd-ac3358da23cb
md"""
And now call the functions for demand and elasticity and store them in the data frame
"""

# ╔═╡ bb612988-5e14-4d2d-ab4d-fe143161d06e
begin # Set the time share = hours/total_hours; the demand function; and elasticity function
	df.time_share .= df.hours ./ sum(df.hours)
	df.D = D((df.a, b, p))
	df.elasticity .= ε_w((b,p,df.D, df.time_share))
end;

# ╔═╡ 72079536-4c38-4c01-97c3-0d6e7f5756fb
begin 
	df.served_demand .= S((df.a, b, p))
	df.VOLL = df.served_demand ./ df.D
end # returns VOLL for the chosen theta, prices, etc

# ╔═╡ 6ba91529-d230-4343-8994-36868e233f9a
# the mean weighted VOLL
mean((df.VOLL .- (df.fuel_cost .+ df.cap_cost)) .* df.time_share)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[compat]
DataFrames = "~1.3.2"
PlutoUI = "~0.7.34"
StatsBase = "~0.33.16"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

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

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "44c37b4636bc54afac5c574d2d02b625349d6582"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.41.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

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

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

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

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

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

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

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

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "0b5cfbb704034b5b4c1869e36634438a047df065"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

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

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

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

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c3d8ba7f3fa0625b062b82853a7d5229cb728b6b"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.1"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8977b17906b0a1cc74ab2e3a05faa16cf08a8291"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.16"

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

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─c160c0bc-0b77-4206-a97c-a0da3d662939
# ╟─ce53326a-2a8a-4bb2-9b5f-c1ec1c2114a7
# ╟─5f6503f1-a50c-454a-baca-55a898146f57
# ╟─88aea7bd-3be8-4a21-8291-8e8d3a3aeca9
# ╟─a28c77ef-e445-4ac2-b183-c4da98f3a4c2
# ╟─ebda0a22-0d15-4eca-b7c8-72a80e27819b
# ╟─d3818f0c-2223-4e04-b5f4-50e632087f61
# ╟─0deaef71-dc51-4aa1-943d-933070f8f89b
# ╟─7b23c101-6b8f-4ecf-878b-0c6b44c3302d
# ╟─ff0d78ba-7597-476c-8549-bfd1577e092c
# ╟─d26c71f8-1a0b-49f8-baed-8dd808dfc08e
# ╟─94da671b-ca21-487e-8b80-d43e82bdbb3e
# ╟─83c7f56e-defe-4f72-9e51-44bdba2b1a0c
# ╟─64bcbe4e-7c16-4959-ad60-d7c64142b4f4
# ╠═d52ae61d-6279-4da2-97dd-934d0f97d4b2
# ╟─b1c42a12-f33c-4c5f-a94d-2867c83b197a
# ╟─a14f7ba2-04be-44ad-853c-9c4732d45109
# ╠═6c2648a8-e894-4f35-be28-d3f8aab665a1
# ╟─a88188b7-6b29-47cf-a34c-71473509cc79
# ╟─a53208d1-7a33-4950-8742-d4430fe621ae
# ╠═2edf0470-88f7-11ec-35a5-c77833b3aad4
# ╟─a0186edd-8afa-46db-8b34-19e275214b5c
# ╠═a7a3a395-a91a-40d2-80da-3e639d3c39b9
# ╟─2dc1da6e-3da5-4681-b7bd-ac3358da23cb
# ╠═bb612988-5e14-4d2d-ab4d-fe143161d06e
# ╠═72079536-4c38-4c01-97c3-0d6e7f5756fb
# ╠═6ba91529-d230-4343-8994-36868e233f9a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
