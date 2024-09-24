
**Doublet_Quantifier**
==================


This code is used as a method to quantify these doublet absorption lines, fitting two absorptions simultaneously, and provide parameters for exploring
aspects of both the shape and size of these absorption lines. The key model in this method used here is the __pseudo-Voigt function__, which is a linear combination of a Gaussian function and Lorentzian profile. 

This work was presented at AAS 243rd meeting as an [iPoster](https://aas242-aas.ipostersessions.com/?s=FA-8F-E2-DC-BE-FF-45-08-3E-48-E5-33-FC-E3-B1-E6).

Work done at Brown Dwarf New York City (BDNYC) and American Museum of Natural History (AMNH)

To cite this work see [CITATION](https://github.com/Ajb2307/Doublet_Quantifier/blob/main/CITATION)


Install Instructions 
===================
First you need to clone this repository, and open the folder.
```
git clone https://github.com/BDNYC/Doublet_Quantifier.git
cd Doublet_Quantifier
```
Next, it needs to be installed.
```
pip install
```



pseudo-Voigt
============
The pseudo-Voigt ($pV$) is a linear combination of the Gaussian distribution ($G$) and a Lorentzian function ($L$) using the mixing coefficient $η$, where $0 \le η \le 1$, where when η=1, the line is a Gaussian distribution and η=0 indicates a Lorentzian distribution.

**Gaussian distribution:** $G(x; σ, μ) = \frac{1}{σ \sqrt{2 π} }e^{-\frac{1}{2} (\frac{x-μ}{σ})^2}$

**Lorentzian function:** $L(x; γ, μ) = \frac{γ}{2 π} \frac{1}{(x-μ)^2 + (γ/2)^2}$

**pseudo-Voigt:** $pV(x; η, A, Γ, μ) = I(η, A, Γ)\left[ η G(x; σ(Γ), μ) + (1-η)L(x; γ(Γ), μ) \right]$

<img src="https://github.com/Ajb2307/Doublet_Quantifier/blob/main/readme-images/pV.png" alt="Pseudo-Voigt with Guassian and Lorentzian plot" width="500">

We use this because it is much simpler than the traditional voigt equation, which is a convolution of a Guassian and Lorentzian disrtibution.

The other three pseudo-Voigt input variables correspond to the maximum depth ($A$), the full-width half-maximum (FWHM, $Γ$), and the center of the distribution ($μ$). These variables are manipulated in the functions $I(η, A, Γ)$ Equation, $σ(Γ)$ Equation, and $γ(Γ)$ when they are inputted into the pseudo-Voigt. 

**Guassian Sigma:** $σ(Γ) = \frac{Γ}{\sqrt{2 \log(2)}}$ 

**Lorenzian Gamma:** $γ(Γ) = 2Γ$ 

**Intensity function:** $I(η, A, Γ) = A \left[ \frac{η}{σ(Γ) \sqrt{2 \pi}} + (1 - η)\frac{2}{\pi γ(Γ)} \right]^{-1}$

<img src="https://github.com/Ajb2307/Doublet_Quantifier/blob/main/readme-images/pV-properties.png" alt="Pseudo-Voigt Properties plotted" width="500">

Fitting the pseudo-Voigt
========================
To fit the doublets, this work uses Python to fit both doublets simultaneously using the pseudo-Voigt model. The code is designed to automate much of the model. The function first requires three inputs from the user: the spectrum as an array, the wavelength range for fitting the continuum, and an approximate wavelength range where the doublet exists.

<img src="https://github.com/Ajb2307/Doublet_Quantifier/blob/main/readme-images/user-inputs.png" alt="vizualized user inputs" width="500">

Then, code does a preliminary fit using a linear continuum fitted to the continuum region. It subtracts this continuum from the spectra before simultaneously fitting two Gaussian distributions to the doublet. The fit is accomplished using SciPy’s curve fit, which optimizes for linear least squares regression. The results of this fit are used to seed a more complex model fit. 

<img src="https://github.com/Ajb2307/Doublet_Quantifier/blob/main/readme-images/preliminary-fit.png" alt="Preliminary Fit plotted" width="500">

Using the preliminary fit, we can constrain the absorption region. Since the absorption region is not stagnate, it can adjust for minor discrepancies in wavelength calibration. The absorption region is determined to be within five standard deviations ($σ$) of the center of the distribution ($μ$), and anything outside of this is considered part of the continuum.

<img src="https://github.com/Ajb2307/Doublet_Quantifier/blob/main/readme-images/updated-bounds.png" alt="Plot of Updated Bounds" width="500">

The next step is fitting a linear continuum to sit on top of the detailed continuum as not to miss the wings of the absorption feature. This is done by binning the continuum region (creating small groups of the wavelength values) and taking the maximum flux within each “bin.” The continuum is then fit to these maximum flux values using SciPy to fit a linear model for the continuum.

<img src="https://github.com/Ajb2307/Doublet_Quantifier/blob/main/readme-images/continuum-fit.png" alt="Continuum Fit plotted" width="500">

The linear continuum is subtracted from the data’s flux, and the pseudo-Voigt fits the data. The starting parameters are chosen utilizing the preliminary fit and data features. The preliminary fit is used as the starting parameter for with the FWHM, $Γ = 2 σ \sqrt{2 \ln{(2)}}$ and the center of each distribution ($μ$) for each feature. The minimum value in the absorption region is used as the starting depth (d). The starting value of the mixing coefficient $η = 0.5$ to initialize the function without bias towards Gaussian distribution or the Lorentzian function.

<img src="https://github.com/Ajb2307/Doublet_Quantifier/blob/main/readme-images/starting-fit.png" alt="Starting Fit plotted" width="500">

The parameters of the fit are optimized twice using SciPy's curve fit, optimizing for least squares regression. The first fit optimizes all the parameters of the pseudo-Voigt, and the second time only adjusts the variable $η$, which shifts the function between a Lorentzian and Gaussian distribution. This is because the mixing parameter ($η$) has little bearing on the overall fit and finely tunes the shape of the absorption, so it is important that it is adjusted on its own once the more influential parameters have been chosen. This results in the final model of the absorption doublet that is returned to the user.

<img src="https://github.com/Ajb2307/Doublet_Quantifier/blob/main/readme-images/final-fit.png" alt="Final Fit plotted" width="500">


Fitted Spectra Examples
=======================
This method has been developed using seven L dwarfs with medium resolution (R~3000) identified in the [SIMPLE Archive](simple-bd-archive.org) and fits well in all but the noisiest of the spectra. This method characterizes the shape and size of these features. It is reliable on a variety of signal-to-noise ratios and is flexible to varying precision of wavelength calibration. 
The method detailed in this paper adds additional numerical representations of the alkali-line doublets, which can be used to further analyze how different parameters in brown dwarfs affect the spectra. However, while this is the intended use of this Python script, nothing limits its use on absorption doublets to any specific case.

<img src="https://github.com/Ajb2307/Doublet_Quantifier/blob/main/readme-images/example-fit1.png" alt="two examples of fitted spectra" width="800">
<img src="https://github.com/Ajb2307/Doublet_Quantifier/blob/main/readme-images/example-fit2.png" alt="two examples of fitted spectra" width="800">

