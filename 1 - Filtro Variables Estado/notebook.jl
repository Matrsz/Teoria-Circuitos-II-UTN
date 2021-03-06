### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 10aba104-a93f-11eb-29ed-4136efb69580
begin
    import Pkg
    Pkg.activate(mktempdir())
    Pkg.add([
        Pkg.PackageSpec(name="Latexify"),
        Pkg.PackageSpec(name="LaTeXStrings"),
        Pkg.PackageSpec(name="Plots"),
		Pkg.PackageSpec(name="ColorSchemes"),
		Pkg.PackageSpec(name="PlotThemes"),
		Pkg.PackageSpec(name="SymEngine"),
		Pkg.PackageSpec(name="PlutoUI")
    ])
    using Latexify
    using LaTeXStrings
	using Plots
	using ColorSchemes
	using PlotThemes
	using SymEngine
	using PlutoUI
end

# ╔═╡ 1a671d44-ac1c-4871-8bcd-7d16642f34e0
begin
	Pkg.add(url="https://github.com/Pocket-titan/DarkMode")
    import DarkMode
end

# ╔═╡ f5b9db94-6ab4-4abb-ac29-6b03ac4dcb16
begin
	md"# Filtro Variables de Estado
	
Tomamos el siguiente circuito, con el objetivo de analizar y simular sus transferencias a V1, V2, V3."
end

# ╔═╡ 1c78f0c7-8485-416f-8789-716571d875b6
md"Planteando las ecuaciones de nodos, se llega al siguiente sistema de 4 ecuaciones con 4 incógnitas y Vi parámetro"

# ╔═╡ 9978620e-7bcc-4c13-b3d9-0dc714e95dae
begin
	nodo1 = :(1/R_4*(V_1-V_X)+1/R_3*(V_3-V_X) = 0)
	nodo2 = :(1/R_1*V_1 + jwC_2*V_3 = 0)
	nodo3 = :(1/R_2*V_2 + jwC_3*V_3 = 0)
	nodo4 = :(1/R_5*(V_i-V_X) + 1/R_6*(V_2-V_X) = 0)
	latexraw(nodo2)
	L"$\begin{aligned}
	\tfrac{1}{R_{4}} \cdot \left( V_{1} - V_{X} \right) + \tfrac{1}{R_{3}} \cdot \left( V_{3} - V_{X} \right) &= 0\\
	\tfrac{1}{R_{1}} \cdot V_{1} + jwC_{2} \cdot V_{2} &= 0\\
	\tfrac{1}{R_{2}} \cdot V_{2} + jwC_{3} \cdot V_{3} &= 0\\
	\tfrac{1}{R_{5}} \cdot \left( V_{i} - V_{X} \right) + \tfrac{1}{R_{6}} \cdot \left( V_{2} - V_{X} \right) &= 0
	\end{aligned}$"
end

# ╔═╡ 7c4cd737-ea86-4d02-92bf-774f8fb04ad0
md"Expresado matricialmente:" 

# ╔═╡ 248a48d7-493b-4506-ae69-3fdbbf82f7c9
L"$\begin{bmatrix}
	\tfrac{1}{R_{4}} + \tfrac{1}{R_3} & -\tfrac{1}{R_4} & 0 & -\tfrac{1}{R_3}\\
	0 & \tfrac{1}{R_1} & jwC_2 & 0\\
	0 & 0 & \tfrac{1}{R_2} & jwC_3\\
	\tfrac{1}{R_{5}} + \tfrac{1}{R_5} & 0 & -\tfrac{1}{R_6} & 0\\
\end{bmatrix}
\begin{bmatrix}
V_X\\ V_2\\ V_3\\ V_4\\
\end{bmatrix}
=
\begin{bmatrix}
0\\ 0\\ 0\\ \tfrac{1}{R_5}\\
\end{bmatrix} V_i$"


# ╔═╡ 44009b12-cbe2-4d18-ba4a-9b4bcf71f465
s = symbols(:s)

# ╔═╡ c5d3db4a-7931-49ae-86dc-a8db73ae2ba1
md"(La funcionalidad interactiva de Pluto Notebooks es todavía incompatible con Github Pages)"

# ╔═╡ fe9697e4-135a-4dd2-9bf1-3d47f6a9da0e
#Rs = @bind R Slider(1:100; default=1, show_value=true)
Rs = 1

# ╔═╡ e01e76cf-6384-4534-8382-eea9a040e25d
#Cs = @bind C Slider(1:100; default=100, show_value=true)
Cs = 100

# ╔═╡ f25d8bcd-5887-4332-b54e-d9fb540825f2
begin
	R1 = Rs*1e3;
	R2 = Rs*1e3;
	R3 = Rs*1e3;
	R4 = Rs*1e3;
	R5 = Rs*1e3;
	R6 = Rs*1e3;
	C2 = Cs*1e-9;
	C3 = Cs*1e-9;
	md"Usando el paquete SymEngine.jl, se puede llegar a una resolución numérica-simbólica, definiendo la variable simbólica s y los valores numéricos de R y C. Los selectores permitirían asignar valores a R en kΩ y C en nF para observar su efecto en la transferencia."
end

# ╔═╡ 71d5083b-ba6c-4a77-a0c9-3b79a282e636
md"En función de eso, se instancia la matriz de admitancias Y y el vector de parámetros i."

# ╔═╡ 1092d595-c1db-4a08-9d18-dacc808e8ae9
Y = [1/R4+1/R3 -1/R4 0    -1/R3;
     0          1/R1 s*C2  0;
     0          0    1/R2  s*C3;
     1/R5+1/R5  0   -1/R6  0]

# ╔═╡ e73c6a27-f6ec-4e05-8fc2-5e277a6c05ce
i = [0; 0; 0; 1/R3]

# ╔═╡ 5a7e24aa-fefb-4e04-8440-2598d3004165
begin
	H = expand.(inv(Y)*i)
	md"SymEngine.jl permite resolver el sistema YH = i semi-simbólicamente, obteniendo así el vector de transferencias H = [Hx; H1; H2; H3] = inv(Y) i, de donde se pueden extraer las transferencias."
end

# ╔═╡ 2324a9f8-a187-4435-81e8-df2c4b44eb90
begin
	H1(s) = H[2]
	H2(s) = H[3]
	H3(s) = H[4]
end

# ╔═╡ 4b87aa13-e337-4512-b477-8e125da89e39
md"Las transferencias corresponden respectivamente a H1 pasa altos, H2 pasa banda, y H3 pasa bajos." 

# ╔═╡ eb0ea7cc-49f9-4e15-9382-193b52ad485c
latexstring("H_1(s) ="*latexraw(H1(s)))

# ╔═╡ 8299fea2-4f51-45b9-a98a-1fb11b511d80
latexstring("H_2(s) ="*latexraw(H2(s)))

# ╔═╡ 1f845ecc-353a-42c8-80dc-159bed6d1652
latexstring("H_3(s) ="*latexraw(H3(s)))

# ╔═╡ 5e678738-638b-44c3-9d6a-df86c3788e01
md"Definiendo la variable simbólica w, SymEngine.jl permite sustituir s = jw en la transferencia H(s) para visualizar las transferencias como diagrama de Bode."

# ╔═╡ eea7c707-8c2f-48e4-86d9-4ede8b111fab
begin
	w = symbols(:w)
	H1w(w) = subs(H1(s), s => w*im)
	H2w(w) = subs(H2(s), s => w*im)
	H3w(w) = subs(H3(s), s => w*im)
end

# ╔═╡ b0fec057-d05b-4a9b-835e-cd267465ff0a
md"A continuación se presentan las gráficas de |H(w)|. Se verifica que H1 corresponde a la salida pasa altos, H2 a la salida pasa banda, y H3 a la salida pasa bajos." 

# ╔═╡ 44bac19d-b875-478e-8a68-2b640765f8cf
begin
	theme(:juno)
	pH1 = plot(abs(H1w(w)), 1, 1e8, xaxis=:log, yaxis=:log, label="H1")
	pH2 = plot(abs(H2w(w)), 1, 1e8, xaxis=:log, yaxis=:log, label="H2")
	pH3 = plot(abs(H3w(w)), 1, 1e8, xaxis=:log, yaxis=:log, label="H3", xlabel="w")
	plot(pH1, pH2, pH3, layout = (3,1))	
end

# ╔═╡ 1df8580c-d92d-4c17-b9f7-d934379c80ad
md"Ahora se procede a una resolución puramente simbólica de la transferencia del filtro, definiendo R y C como variables. Por simplicidad del resultado final suponemos resistencias y capacitores de igual valor." 

# ╔═╡ cfee3c91-9ec7-4e1f-9fda-f7b9fb2407e0
@vars R C

# ╔═╡ 51bee854-294a-44c3-89ae-392969e53c08
Ys = [1/R+1/R -1/R  0   -1/R;
      0        1/R  s*C  0;
      0        0    1/R  s*C;
     1 /R+1/R  0   -1/R  0]

# ╔═╡ 690f4c2c-a454-4c18-b798-13f4e2c08eb8
is = [0; 0; 0; 1/R]

# ╔═╡ 8d914ac6-b6fd-4ab9-a7d4-40645205c3ea
begin
	Hs = expand.(inv(Ys)*is)
	Hs1(s) = Hs[2]
	Hs2(s) = Hs[3]
	Hs3(s) = Hs[4]
	md"Con un procedimiento idéntico al anterior se puede obtener una expresión puramente simbólica de las transferencias H1, H2, H3:"
end

# ╔═╡ 0d8fa251-0e1e-4f24-9e4b-5eab54cb483d
latexstring("H_1(s) ="*latexraw(Hs1(s)))

# ╔═╡ 3ecde8b1-5b0c-4a78-8768-46a5cdf25dff
latexstring("H_2(s) ="*latexraw(Hs2(s)))

# ╔═╡ 28fe2ae9-9f71-4dfb-a8a0-bcb925498978
latexstring("H_3(s) ="*latexraw(Hs3(s)))

# ╔═╡ 35ee3ea7-465e-4890-8bd9-98066ddfb10a
md"O bien, escrito de la forma habitual:"

# ╔═╡ 7642f158-b9ab-4ca5-af98-0ca8fb9c053c
	L"$\begin{aligned}
	H_1(s) &= \frac{R^2C^2s^2}{R^2C^2s+RCs+1}\\[1em]
	H_2(s) &= \frac{RCs}{R^2C^2s+RCs+1}\\[1em]
	H_1(s) &= \frac{1}{R^2C^2s+RCs+1}\\
	\end{aligned}$"

# ╔═╡ 8a097b97-65bb-4b9f-87c2-25226ab2a7fc
md"Es decir, resultan tres filtros con frecuencia de corte wc=RC, que en el caso particular correspondiente a la resolución numérica de R = 1kΩ y C = 100nF retorna una frecuencia wc = 10000 rad/s, consistente con lo visto en los diagramas de Bode."

# ╔═╡ 6ec8a49a-fb90-493a-8301-9aa23c482c25
begin
	struct Wow
		filename
	end

	function Base.show(io::IO, ::MIME"image/png", w::Wow)
		write(io, read(w.filename))
	end
end

# ╔═╡ 98c4b949-7464-4690-92cc-93e205a35072
Wow("Images/circuito.png")

# ╔═╡ e8fca5e5-2263-4e73-9481-292a528961c1
DarkMode.enable(theme="nord")

# ╔═╡ Cell order:
# ╟─f5b9db94-6ab4-4abb-ac29-6b03ac4dcb16
# ╟─98c4b949-7464-4690-92cc-93e205a35072
# ╟─1c78f0c7-8485-416f-8789-716571d875b6
# ╟─9978620e-7bcc-4c13-b3d9-0dc714e95dae
# ╟─7c4cd737-ea86-4d02-92bf-774f8fb04ad0
# ╟─248a48d7-493b-4506-ae69-3fdbbf82f7c9
# ╟─f25d8bcd-5887-4332-b54e-d9fb540825f2
# ╠═44009b12-cbe2-4d18-ba4a-9b4bcf71f465
# ╟─c5d3db4a-7931-49ae-86dc-a8db73ae2ba1
# ╠═fe9697e4-135a-4dd2-9bf1-3d47f6a9da0e
# ╠═e01e76cf-6384-4534-8382-eea9a040e25d
# ╟─71d5083b-ba6c-4a77-a0c9-3b79a282e636
# ╠═1092d595-c1db-4a08-9d18-dacc808e8ae9
# ╠═e73c6a27-f6ec-4e05-8fc2-5e277a6c05ce
# ╟─5a7e24aa-fefb-4e04-8440-2598d3004165
# ╠═2324a9f8-a187-4435-81e8-df2c4b44eb90
# ╟─4b87aa13-e337-4512-b477-8e125da89e39
# ╠═eb0ea7cc-49f9-4e15-9382-193b52ad485c
# ╠═8299fea2-4f51-45b9-a98a-1fb11b511d80
# ╠═1f845ecc-353a-42c8-80dc-159bed6d1652
# ╟─5e678738-638b-44c3-9d6a-df86c3788e01
# ╠═eea7c707-8c2f-48e4-86d9-4ede8b111fab
# ╟─b0fec057-d05b-4a9b-835e-cd267465ff0a
# ╟─44bac19d-b875-478e-8a68-2b640765f8cf
# ╟─1df8580c-d92d-4c17-b9f7-d934379c80ad
# ╟─cfee3c91-9ec7-4e1f-9fda-f7b9fb2407e0
# ╟─51bee854-294a-44c3-89ae-392969e53c08
# ╟─690f4c2c-a454-4c18-b798-13f4e2c08eb8
# ╟─8d914ac6-b6fd-4ab9-a7d4-40645205c3ea
# ╟─0d8fa251-0e1e-4f24-9e4b-5eab54cb483d
# ╟─3ecde8b1-5b0c-4a78-8768-46a5cdf25dff
# ╟─28fe2ae9-9f71-4dfb-a8a0-bcb925498978
# ╟─35ee3ea7-465e-4890-8bd9-98066ddfb10a
# ╟─7642f158-b9ab-4ca5-af98-0ca8fb9c053c
# ╟─8a097b97-65bb-4b9f-87c2-25226ab2a7fc
# ╟─6ec8a49a-fb90-493a-8301-9aa23c482c25
# ╟─10aba104-a93f-11eb-29ed-4136efb69580
# ╟─1a671d44-ac1c-4871-8bcd-7d16642f34e0
# ╟─e8fca5e5-2263-4e73-9481-292a528961c1
