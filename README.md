Doublet_Quantifier
==================

This code is used as a method to quantify these doublet absorption lines, fitting two absorptions simultaneously, and provide parameters for exploring
aspects of both the shape and size of these absorption lines. The key model in this method used here is the __pseudo-Voigt function__, which is a linear combination of a Gaussian function and Lorentzian profile. 

This work was presented at AAS 243rd meeting as an [iPoster](https://aas242-aas.ipostersessions.com/?s=FA-8F-E2-DC-BE-FF-45-08-3E-48-E5-33-FC-E3-B1-E6).

Work done at Brown Dwarf New York City (BDNYC) and American Museum of Natural History (AMNH)

pseudo-Voigt
============
The pseudo-Voigt ($pV$) is a linear combination of the Gaussian distribution ($G$) and a Lorentzian function ($L$) using the mixing coefficient $η$, where $0 \le η \le 1$, where when η=1, the line is a Gaussian distribution and η=0 indicates a Lorentzian distribution.

**Gaussian distribution:** $G(x; σ, μ) = \frac{1}{σ \sqrt{2 π} }e^{-\frac{1}{2} (\frac{x-μ}{σ})^2}$

**Lorentzian function:** $L(x; γ, μ) = \frac{γ}{2 π} \frac{1}{(x-μ)^2 + (γ/2)^2}$

**pseudo-Voigt:** $pV(x; η, A, Γ, μ) = I(η, A, Γ)\left[ η G(x; σ(Γ), μ) + (1-η)L(x; γ(Γ), μ) \right]$

We use this because it is much simpler than the traditional voigt equation, which is a convolution of a Guassian and Lorentzian disrtibution.

The other three pseudo-Voigt input variables correspond to the maximum depth ($A$), the full-width half-maximum (FWHM, $Γ$), and the center of the distribution ($μ$). These variables are manipulated in the functions $I(η, A, Γ)$ Equation, $σ(Γ)$ Equation, and $γ(Γ)$ when they are inputted into the pseudo-Voigt. 

$σ(Γ) = \frac{Γ}{\sqrt{2 \log(2)}}$ 

$γ(Γ) = 2Γ$ 

$I(η, A, Γ) = A \left[ \frac{η}{σ(Γ) \sqrt{2 \pi}} + (1 - η)\frac{2}{\pi γ(Γ)} \right]^{-1}$




