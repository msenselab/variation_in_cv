

_This is the repository for the study "Variation in the 'coefficient of variation': rethinking the violation of the scalar proerpty in time-duration judgments", Yue Ren, Fredrik Allenmark, Hermann J. Müller, Zhuanghua Shi_

1. Please download all the files from the this repository in one local folder, and designate the same folder as the current working directory. 

2. Running the file 'Main_Figure.R' independently returns:  1), all the figure plots (Figure 2 and 3) in the manuscript, and 2), data source including raw data, basic statistics (e.g., grand-mean), and modeling fittings. 

(tips: Please collapse all the chunks (separate with ----#) for a clearer structure.)

3. The whole procedure of the model fitting can be found in 'Bayesian_model_final.R', which can also be run independently. 

The model assumes that Bayesian integration happens on a logarithmic scale representation
of the durations, followed by a transformation back to linear scale for the reproduction
during which some additional decision and motor noise is introduced. The Bayesian integration
on logarithmic scale on its own would predict constant CV, but the introduction of additional
duration-independent noise after transforming back to linear scale results in the predicted
CV being larger for the shortest durations.

For the blocked condition there is additionally a model comparison between 32 different
models which differ in terms of whether sig_p, sig_s, sig_m and res can differ between short, 
medium and long duration blocks as well as in whether the bias "res" is added as a shift of 
the mean of the prior, or later as a bias of the reproduction