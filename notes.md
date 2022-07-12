# Deconvolution of Enhancers in Mixed Nascent Sequencing Samples

This is the start of this project's lab notebook.

Basic outline:
- [B] TFit or dREG on all the constituent samples
- [B] muMerge keeping singletons to get them all together
- [C] DESeq2 likelihood-ratio to run automagically between samples and infer differences
- [A] Run nu-SVM to infer parameters
  - [B] Use l-curve optimization on the SVM
