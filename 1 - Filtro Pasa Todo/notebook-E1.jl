### A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

# ╔═╡ 53d9d2c0-0369-403e-8378-858663608a35
begin
    import Pkg
    Pkg.activate(mktempdir())
    Pkg.add([
		Pkg.PackageSpec(name="Revise"),
        Pkg.PackageSpec(name="Latexify"),
        Pkg.PackageSpec(name="LaTeXStrings"),
        Pkg.PackageSpec(name="ControlSystems"),
		Pkg.PackageSpec(name="ColorSchemes"),
		Pkg.PackageSpec(name="PlotThemes")
    ])
	using Revise
    using Latexify
    using LaTeXStrings
	using ControlSystems
	using ColorSchemes
	using PlotThemes
end

# ╔═╡ 451411e9-6eb1-4ad1-97cd-946300511c0d
begin
	Pkg.add(url="https://github.com/Pocket-titan/DarkMode")
    import DarkMode
end

# ╔═╡ 7db3794b-9194-428c-a745-8ec848acadb5
begin
	using Plots
	theme(:juno)
end

# ╔═╡ 550a6097-7c37-41dc-aa8c-0da7b523f31a
begin
	md"# Teoría de Circuitos II - Filtro Pasa Todo
	
Tomamos el siguiente circuito, con el objetivo de analizar y simular su transferencia."
end

# ╔═╡ 67c80a93-838d-4f72-b917-75ff417e1d71
md"Si se usa el modelo OpAmp ideal $:\, V^+=V^-=V_x, \quad Z_{in}\rightarrow\infty$ se encuentra la siguiente ecuación matricial aplicando ley de nodos a $V^+$y $V^-$:"

# ╔═╡ 70bb620a-b3d3-4d48-a591-f2e3abbb5f59
begin
	ex1 = :([1/R_1+1/R_2 -1/R_2; 1/R_1+jwC 0])
	ex2 = :([V_x; V_2])
	ex3 = :([1/R_1; jwC])
	latexstring(latexraw(ex1),latexraw(ex2),"=",latexraw(ex3),"V_1")
end

# ╔═╡ b9353662-c3ca-460d-8724-ff01d8692b32
md"La transferencia $H(jw) = \frac{V_2}{V_1}$ se puede obtener despejando la ecuación y resolviendo para $V_2$."

# ╔═╡ ec139584-2a0a-4583-8392-ee5bcd7096f7
begin
	ex4 = :([1/R_1+1/R_2 -1/R_2; 1/R_1+jwC 0]^(-1))
	latexstring(latexraw(ex2),"=",latexraw(ex4),latexraw(ex3),"V_1")
end

# ╔═╡ 14dfa10b-ab0d-4b6b-b3b9-aad4838c157b
md"$\vdots$ Contemplando $R_1 = R_2$:"

# ╔═╡ 8eaa07ab-f6f1-45aa-a929-f25958bfef73
begin
	ex5 = :(H(jw)=(jCw-1/(R_3))/(jCw+1/(R_3)))
	latexify(ex5)
end

# ╔═╡ 67d13c95-9f5d-4753-90ed-b15612b89945
md"La simulación numérica la vamos a realizar en el lenguaje Julia, usando el módulo ControlSystems.jl"

# ╔═╡ 97467fe3-f873-40a0-a402-f4c2d398b956
begin
	    R1 = 1e3;
	    R2 = 1e3;
	    C  = 1e-6;
	    R3 = 1e3;
	md"Una vez definidas las variables R1, R2, C, R3, definimos los coeficientes del numerador y denominador."
end

# ╔═╡ 3f47cbfb-ccb5-4f19-99fb-387ccdc26554
num = [C/R2, -1/(R1*R3)]

# ╔═╡ 536d941c-24e3-454d-b8af-1da615cee757
den = [C/R2,  1/(R2*R3)]

# ╔═╡ 578b5413-6af8-429b-8e6f-4581704a21e6
md"ControlSystems.tf retorna H como un objeto TransferFunction:"

# ╔═╡ 5d68641a-ac0b-48c8-b7aa-9f4798f1f729
H = tf(num, den)

# ╔═╡ 0b356ae9-a7a0-4e19-add1-1d7498c10869
md"ControlSystems.bodeplot toma el obtjeto TransferFunction y obtiene la respuesta en frecuencia,que efectivamente corresponde a un filtro pasa todo." 

# ╔═╡ 006447ea-7984-4eb2-916d-e85181b639c3
bodeplot(H, palette=:seaborn_muted6)

# ╔═╡ 77224aa9-e2cc-493c-a810-4affe4b33428
md"ControlSystems.rlocusplot toma el objeto TransferFunction y obtiene el diagrama de polos y ceros, que evidencia la estabilidad del sistema y hace visible la condición de filtro pasa todo: simetría de polos y ceros respecto al eje jw."

# ╔═╡ 2162cb78-efdf-4b65-bed0-e87b2b28d829
rlocusplot(H; aspect_ratio=:equal, framestyle=:origin, palette=:seaborn_muted6)

# ╔═╡ cd8440fe-d58d-4c2f-a233-2f1aaec2c647
md"Finalmente, volviendo al modelo en LTSpice comparamos los resultados de la simulación numérica con una simulación circuital"

# ╔═╡ c5d02935-1c7c-47bd-8b52-d0f238c8cd9b
md"La simulación de LTSpice, a diferencia del modelo ideal, implementa el modelo integrador del amplificador operacional. El resultado es el siguiente:" 

# ╔═╡ 26701700-7abf-40a6-a11f-776684f48ba7
md"Es un comportamiento equivalente a la simulación numérica a bajas frecuencias, filtro pasa todo, pero el comportamiento difiere a altas frecuencias donde el efecto del capacitor del modelo integrador es apreciable, y la transferencia real disminuye." 

# ╔═╡ b69f9cbc-b116-4206-9ede-2a3a2d6bc135
begin
	struct Wow
		filename
	end

	function Base.show(io::IO, ::MIME"image/png", w::Wow)
		write(io, read(w.filename))
	end
end

# ╔═╡ 2d83491c-b238-40d2-903a-e6368089bdd3
Wow("Images/circuito.png")

# ╔═╡ ecd5ba2f-0337-4e8d-bdb0-4e8e2c15607b
Wow("Images/circuito.png")

# ╔═╡ 0193b0ca-26ff-45ce-bdaa-650a5fc807dd
Wow("Images/sim.png")

# ╔═╡ 40ed78f5-ab4a-4fcc-b204-a1cec7f8595a
DarkMode.enable(theme="nord")

# ╔═╡ Cell order:
# ╟─550a6097-7c37-41dc-aa8c-0da7b523f31a
# ╟─2d83491c-b238-40d2-903a-e6368089bdd3
# ╟─67c80a93-838d-4f72-b917-75ff417e1d71
# ╟─70bb620a-b3d3-4d48-a591-f2e3abbb5f59
# ╟─b9353662-c3ca-460d-8724-ff01d8692b32
# ╟─ec139584-2a0a-4583-8392-ee5bcd7096f7
# ╟─14dfa10b-ab0d-4b6b-b3b9-aad4838c157b
# ╟─8eaa07ab-f6f1-45aa-a929-f25958bfef73
# ╟─67d13c95-9f5d-4753-90ed-b15612b89945
# ╟─97467fe3-f873-40a0-a402-f4c2d398b956
# ╠═3f47cbfb-ccb5-4f19-99fb-387ccdc26554
# ╠═536d941c-24e3-454d-b8af-1da615cee757
# ╟─578b5413-6af8-429b-8e6f-4581704a21e6
# ╠═5d68641a-ac0b-48c8-b7aa-9f4798f1f729
# ╟─0b356ae9-a7a0-4e19-add1-1d7498c10869
# ╠═006447ea-7984-4eb2-916d-e85181b639c3
# ╟─77224aa9-e2cc-493c-a810-4affe4b33428
# ╠═2162cb78-efdf-4b65-bed0-e87b2b28d829
# ╟─cd8440fe-d58d-4c2f-a233-2f1aaec2c647
# ╟─ecd5ba2f-0337-4e8d-bdb0-4e8e2c15607b
# ╟─c5d02935-1c7c-47bd-8b52-d0f238c8cd9b
# ╟─0193b0ca-26ff-45ce-bdaa-650a5fc807dd
# ╟─26701700-7abf-40a6-a11f-776684f48ba7
# ╟─b69f9cbc-b116-4206-9ede-2a3a2d6bc135
# ╟─53d9d2c0-0369-403e-8378-858663608a35
# ╟─451411e9-6eb1-4ad1-97cd-946300511c0d
# ╟─40ed78f5-ab4a-4fcc-b204-a1cec7f8595a
# ╟─7db3794b-9194-428c-a745-8ec848acadb5
