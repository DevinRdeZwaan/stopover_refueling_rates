#####################################################################################
# R Code for:

# Mass gain and stopover dynamics among migrating songbirds are linked to seasonal, environmental, and life-history effects

# Authors: DR de Zwaan, A Huang, Q McCallum, K Owen, M Lamont, and W Easton (2022)

#####################################################################################

### Project goals:
During migration, birds must stopover at refueling sites to replenish energy stores, with the rate of refueling linked to stopover length, migration speed, and the timing of critical life-stages, such as breeding. Under optimal migration theory, birds are expected to maximize fuel intake and minimize stopover length (‘time minimization hypothesis’). Using long-term banding data (10 years) from the Iona Island Bird Observatory (IIBO) in southwestern British Columbia, Canada, we assessed spring and fall stopover behaviors among five warbler (Parulidae) and five sparrow species (Passerellidae) that differ in their migration strategies, diet, and foraging niche. Songbirds were captured in both spring and fall, providing a powerful comparative framework in which to test refueling rate and stopover dynamic associations among species and seasons. Specifically, we assessed: (1) species-specific refueling rates between seasons, diet guilds (insectivores, omnivores), and migration distances (short, long), (2) season-specific influences of temperature, precipitation, and habitat productivity on refueling rate, (3) density-dependent effects on refueling performance, and (4) the influence of body mass on departure probability and stopover duration. 


# Table of Contents:

# 1) '01_de Zwaan et al_fuel index analysis'

### Description
Conducts the fuel index analysis and the association between mass gain, departure probability, and stopover duration using birds that have been captured multiple times within a season. This code is specifically for spring migration data, but fall migration analysis can be run on the same code by substituting the spring for fall data.

# 2) '02_de Zwaan et al_mark-recapture code for stopover duration'

### Description
This code: A) calculates lean body mass (LBM) and fuel index from raw banding data, B) estimates the 'true' stopover duration of all birds using capture-recapture data, and C) assesses how fuel index at first capture (a proxy for arrival condition) is associated with length of stay. Again, note that this code demonstrates the work flow for the spring migration analysis only. To conduct the fall migration analysis, simply substitute the spring data with the fall data.
