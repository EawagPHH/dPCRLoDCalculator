# dPCRLoDCalculator

This Rshiny app allows you to easily calculate LoD using two models. Simply upload your QX600 data .csv with the "Sample Description 1" column filled with estimated expected concentrations of your individual measurements.

# **Understanding Limit of Detection (LoD)**

The **Limit of Detection (LoD)** is the lowest concentration of DNA that can be reliably detected, not just occasionally, but with high confidence. In ddPCR, we define this as the concentration where there is a **95% chance** of detecting the target.

In simple terms:

> If we tested many samples at this concentration, we would detect the DNA 95 out of 100 times.

---

# **Why 95%?**

The 95% standard is widely accepted in clinical diagnostics and environmental monitoring. It ensures the false-negative rate stays below 5%.

This standard is recommended in:

**CLSI EP17-A2**  
*Evaluation of Detection Capability for Clinical Laboratory Measurement Procedures*  
https://clsi.org/standards/products/method-evaluation/documents/ep17/

---

# **Two Statistical Models Used to Estimate LoD**

Detection is binary:

- 1 = detected  
- 0 = not detected  

We model the **probability of detection** as concentration increases.

---

## 1️⃣ **Model 1:  Probit Regression**

The Probit model assumes detection probability follows a smooth S-shaped curve.

The mathematical model is:

$$
\Phi^{-1}(p) = \beta_0 + \beta_1 X
$$

What this means in plain English:

- $X$ = concentration
- $p$ = probability of detecting DNA
- $\beta_0$ shifts the curve left or right
- $\beta_1$ controls how steep the curve is

If the slope ($\beta_1$) is steep:
> A tiny increase in DNA concentration causes a large jump in detection probability.

To find the concentration where detection = 95%, we solve:

$$
LoD_{95} = \frac{z_{0.95} - \beta_0}{\beta_1}
$$

You don’t need to compute this manually since the app does it for you.

Conceptually:

> We’re finding the concentration where the model predicts 95% detection.

Reference:  
Forootan et al. 2017  
https://doi.org/10.1016/j.bdq.2017.04.001

---

## 2️⃣ **Model 2: Weibull Model**

Sometimes detection probability does not increase symmetrically.

At very low concentrations, detection may “trail off” more gradually. It's not a perfect S shape.

The Weibull model captures this behavior:

$$
p(X) = 1 - \exp \left[-\left(\frac{X}{\lambda}\right)^k \right]
$$

Plain explanation:

- $\lambda$ controls scale (where the curve sits)
- $k$ controls shape (how curved or stretched it is)

To find LoD:

$$
LoD_{95} = \lambda \left[-\ln(1 - 0.95)\right]^{1/k}
$$

In simple terms:

> We plug 95% into the curve and solve for concentration.

Reference:  
Klymus et al. 2019  
https://doi.org/10.1002/edn3.29

---

## **How Are Model Parameters Estimated?**

Both models use **Maximum Likelihood Estimation (MLE)**.

The likelihood equation is:

$$
L = \prod_{i=1}^{n} p_i^{y_i} (1 - p_i)^{1 - y_i}
$$

You don’t need to follow the math.

What it means:

> The algorithm finds the curve that makes your actual observed results (detected vs not detected) most probable.

It is not guessing but rather it is mathematically optimizing the best fit.

---

## **Why Confidence Intervals Matter**

The LoD is an estimate, not an exact value.

We calculate a 95% confidence interval to show uncertainty.

Interpretation:

> If we repeated this entire experiment many times, 95% of calculated LoDs would fall within this range.

---

## **Probit Confidence Interval (Delta Method)**

The LoD is calculated from model coefficients:

$$
LoD_{95} = \frac{z_{0.95} - \beta_0}{\beta_1}
$$

Because this is a formula involving estimated parameters, we approximate its uncertainty using the **Delta Method**.

In simple terms:

> We estimate how much the error in the slope and intercept affects the final LoD.

---

## **Weibull Confidence Interval (Bootstrapping)**

Instead of using formulas, we use simulation:

1. Randomly resample wells (with replacement)
2. Recalculate LoD
3. Repeat 100 times
4. Take the middle 95% of results

Plain explanation:

> We ask 100 slightly different versions of your dataset what the LoD is.  
> If they mostly agree, confidence is high.

---

## **Model Selection (RSS)**

We compare models using:

$$
RSS = \sum (Observed - Predicted)^2
$$

Plain explanation:

> RSS measures how far the model’s predictions are from your real data.

Lower RSS = better fit.

If the difference between models is very small (< 0.05), we consider them practically equivalent.

---

# **What If One Model Is Recommended Over the Other?**

Sometimes the app will suggest either the Probit or Weibull model based on a lower RSS value.

## **What does that mean?**

It simply means:

> One mathematical curve fits your observed detection data slightly better.

This does **not** mean the other model is “wrong.”  
It only indicates which curve more closely matches the pattern in your data.

---

## **If Probit Is Recommended**

This usually means:

- Detection probability increases smoothly and symmetrically.
- The assay behaves like a typical dose–response curve.
- There is no unusual trailing at very low concentrations.

In this case, you can confidently report the Probit-based LoD.

Probit is widely used in clinical diagnostics and aligns well with traditional regulatory expectations.

---

## **If Weibull Is Recommended**

This usually means:

- Detection probability rises asymmetrically.
- There is a gradual “tail” at very low concentrations.
- The assay may show variability near the detection limit.

In environmental DNA (eDNA) and low-template samples, this behavior is common.

In these situations, the Weibull model may provide a more realistic estimate of LoD.

---

## **What Should You Do?**

If one model has a clearly lower RSS:

✔ Report the recommended model  
✔ Include the 95% confidence interval  
✔ Document which model was used  

If the RSS values are nearly identical:

✔ Either model is acceptable  
✔ Report the simpler model (typically Probit) for transparency  
✔ Note that both models produced similar results  

---

## **Regulatory or Publication Considerations**

For clinical validation:

- Probit is often preferred due to historical precedent.
- However, documenting model comparison strengthens your justification.

For environmental or research applications:

- The Weibull model may better reflect biological detection patterns.

---

## **Best Practice Recommendation**

Always:

1. Report the selected model.
2. Report the LoD with its 95% confidence interval.
3. Document how model selection was performed (RSS comparison).

Transparency is more important than model choice.

The goal is not to choose the “most complex” model —  
it is to choose the one that best represents your data.

---

# **Other Data Quality Requirements**

**Dilutions and detects**  
In order for LoD to be properly calculated, it is recommended to have 7 two-fold dilutions of which **at least one concentration** has a detect rate of 100%, and 2-6 have detect rates below 100%. 

**Minimum 10,000 droplets**  
Below this threshold, the Poisson assumptions behind ddPCR becomes less reliable but it is still usable if you wish, just acknowledge that replicates become less precise.

**Reaction Volume**  
Enter only the volume actually loaded into the QX600 droplet reader (almost always 20).  
Do not include extra pipetting volume to account for pipetting losses.

---
