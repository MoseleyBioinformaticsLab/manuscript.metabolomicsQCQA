
## Analyses that maybe could be addressed by RCHAMetab??

* AN000112:
  * No common sample IDs between SSF and measurements due to the extra stuff tacked onto the measurements sample IDs.
* AN003604;
  * Something goes wrong with the cleaning sample_ids, I think. Yeah, there are "-" and "+" that get converted to just a single "_"
* AN004188:
  * has a "+" at the end, which doesn't get converted properly.
* AN000144:
  * One sample, S00011727, has duplicated measurements.
* AN000139:
  * Has QC samples in the measurements, that are not mentioned in the SSF at all.
* AN000627:
  * Missing the MW block.
* AN003915:
  * Has empty keys in the SSF block. I drop the whole thing because of the empty keys.

## Analyses that we want to highlight and find by some methods in QC/QA
* AN000146:
  * The SSF are incomplete for a large number of samples, and it looks like some of the SSF data should actually be part of the extra data.
* AN000072:
  * Has a bunch of data points that look like they are zero?? But didn't get set to NA when the data was being normalized. Nope, they are all equal to "1.000". Wow. But there are also NA values, so it's very, very weird.
* AN000065:
  * Has 2 healthy, and 19 injured cases.
* AN000123:
  * Had "quality check* for the name of some samples". - this we are catching now.
* AN000070:
  * Study says 4 groups with 75 subjects, but AN000070 only has 3 groups with 30 subjects.
* AN001723:
  * Has negative values. Ugh. I wonder how many things in our data have negative values and we just didn't catch it.
* AN002289:
  * Pretty sure this has a sample with all NA in it. Check if it's real!
* AN003701
  * Might be the same deal with a sample with all NA.
* AN000163:
  * Has no measurements in the mwtab, and no ancillary data.
* AN004749:
  * Has all same values, or missing.
* AN005958:
  * Has ".E" instead of just "E" to denote the exponent, which had to be corrected.