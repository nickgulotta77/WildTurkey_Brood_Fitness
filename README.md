# Fitness consequences of risk-taking by brooding female wild turkeys in the southeastern United States
### Kevin J. Kinser-Mcbee1,†1, Nick A. Gulotta1,†1, Nicholas W. Bakner2, Erin E. Ulrey3, Bradley S. Cohen4, Bret A. Collier5, Michael J. Chamberlain1
### 

####  † Kevin J. Kinser-Mcbee and Nick A. Gulotta are joint first authors
####  1 Warnell School of Forestry and Natural Resources, University of Georgia, Athens, Georgia, USA
####  2 Department of Entomology and Wildlife Ecology, University of Delaware, Newark, DE, USA.
####  3 School of Natural Resources, University of Nebraska-Lincoln, Lincoln, USA
####  3 College of Arts and Sciences, Tennessee Technological University, Cookeville, Tennessee, USA
####  4 Department of Fisheries and Wildlife, Michigan State University, East Lansing, MI 48824



## Abstract
Survival-reproduction trade-offs are well studied in altricial birds, where offspring remain at nests and parental care can be shared between adults. Less is known about these trade-offs in species with uniparental care, where a single parent must navigate heterogeneous landscapes while caring for precocial offspring. In such systems, habitats that benefit offspring may expose parents to greater mortality risk, meaning habitat management focused on reproduction alone could overlook costs to adult survival. In these systems, risk-taking may reflect resource allocation, where risky areas benefit offspring at the cost of adult survival; resource acquisition, where risky areas support both adult and offspring survival; or maladaptive risk-taking, where risky areas reduce both fitness components without compensating benefits. We used GPS movement data from 207 brooding female eastern wild turkeys (Meleagris gallopavo silvestris) monitored across seven study sites in the southeastern United States from 2014-2025 to quantify individual variation in habitat associated risk-taking and evaluate consequences for female and brood survival. We found that the fitness consequences of habitat-associated risk-taking varied across landcover types. Females closer to open landcover had lower survival (β = 0.41, 95% CrI = 0.05, 0.83) but higher brood survival (β = -0.27, 95% CrI = -0.46, -0.09), indicating resource allocation. Females closer to hardwood landcover had higher individual survival and higher brood survival, consistent with resource acquisition. By contrast, females closer to pine landcover had lower individual survival and lower brood survival, suggesting maladaptive risk-taking. Together, these findings reveal an important consideration for habitat management, habitats associated with greater reproductive success may not be those associated with greater adult survival. By linking among-individual differences in habitat-associated risk-taking to multiple fitness components, our framework offers managers of species of conservation concern a means of identifying hidden demographic trade-offs and prioritizing habitats associated with complementary benefits to both adults and offspring.
## General Project File Structure

```
├── README.md                                                          <- The top-level README including general project descriptions
|
├── Behavioral data
│   ├──daily_landcover_brood_2014-2025.csv                              <- Behavioral data ready for modeling
|
├── Survival data
│   ├── brood_survival_ref_2014_2025_FINAL.csv                         <- Survival data ready for modeling
|
├── Rds files
│   ├── L.zip                                                          <- Fitted model file for hardwood landcover
│   ├── L.zip                                                          <- Fitted model file for pine landcover
│   ├── L.zip                                                          <- Fitted model file for open landcover
|
├── R script files
│   ├── ACF_plots.r                                                    <- script to produce ACF plots
│   ├── Univariate_hardwood_pine_open_fitness_.r                       <- Script to produce main results and figures
