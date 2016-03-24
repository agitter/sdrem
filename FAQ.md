# SDREM Frequently Asked Questions

### Contents
- [TF activity timing](#tftiming)
- [Randomization error](#randerror)

## <a id='tftiming'>TF activity timing</a>
**Q:** How do I identify the time points at which TFs are first active?

**A:** Visually, the Key Transcription Factors interface in the DREM display can be used to show TFs the first time they are active on a DREM path based on the SDREM activity score.  This information is not directly saved to a text file, but it can be reconstructed from the ```N.model.activitiesDynamic``` file that SDREM produces.

```N.model.activitiesDynamic```, where *N* is the final iteration of SDREM, is a tab-delimited text file that contains three columns:
- TF identifier
- TF activity score at a particular split node in the DREM paths
- Tree depth of that split node

The tree depth of the split node can be used to identify the first time a TF is active.  The depth is 0-based so if you have time points of 0, 5, 15, 30, 60 minutes the corresponding tree depths are 0, 1, 2, 3, and 4.  By parsing or manually examining ```N.model.activitiesDynamic```, you can identify the first time (minimum tree depth) at which a TF exceeds the activity score threshold you set in the Key Transcription Factors slider.  Note that for activity scores, the slider is used to threshold TFs that have a score >= 10^X.

To illustrate, suppose you set X = 4.5 in the slider.  Then a TF will be shown the first time its activity exceeds 10^4.5 = 31622 on a DREM path.  You can use ```N.model.activitiesDynamic``` to filter out the rows where the activity score is < 31622.  This will give you a list of TF identifiers and tree depths of the split nodes where those TFs attained the high activity score.  You can sort by TF identifiers to determine the smallest tree depth with a high activity score for each TF.  This gives you an index, which can be used to recover the time point.  Using the example time points from above, if a TF exceeds the activity score threshold at split nodes with depth 2, 3, and 3, you would take 2 as the minimum depth (time point index), which corresponds to the 15 minute time point.

The ```N.model.activitiesDynamic``` file can also be used to identify all time points at which a TF is active.  In the example above, the TF is active at split nodes with depth 2 and 3 meaning it is active at 15 and 30 minutes.

## <a id='randerror'>Randomization error</a>
**Q:** Why does SDREM throw a randomization error in ```alg.DREMInterface.randomizeBindingPriors```?

**A:** This can occur when the TF-gene interactions file has multiple rows for the same gene.  SDREM assumes that each TF only appears in one column and each gene only appears in one row.
