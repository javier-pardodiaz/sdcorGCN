from __future__ import print_function
import os
import sys

##########################################################
## This script combines all the STRING's channels subscores
## into the final combined STRING score.
## It uses unpacked protein.links.full.xx.txt.gz as input
## which can be downloaded from the download subpage:
##      https://string-db.org/cgi/download.pl
##########################################################
 
input_file = "4932.protein.links.full.v11.0.txt"

if not os.path.exists(input_file):
    sys.exit("Can't locate input file %s" % input_file)

prior = 0.041

def compute_prior_away(score, prior):

    if score < prior: score = prior
    score_no_prior = (score - prior) / (1 - prior)

    return score_no_prior

header = True
file = open("Yeast_coexpression.txt","w") 
 

for line in open(input_file):

    if header:
        header = False
        continue
    
    l = line.split()
    
    ## load the line
        
    (protein1, protein2,
     neighborhood, neighborhood_transferred,
     fusion, cooccurrence,
     homology,
     coexpression, coexpression_transferred,
     experiments, experiments_transferred,
     database, database_transferred,
     textmining, textmining_transferred,
     initial_combined) = l


    ## divide by 1000

    neighborhood = float(neighborhood) / 1000
    neighborhood_transferred = float(neighborhood_transferred) / 1000
    fusion = float(fusion) / 1000
    cooccurrence =  float(cooccurrence) / 1000
    homology = float(homology) / 1000
    coexpression = float(coexpression) / 1000
    coexpression_transferred = float(coexpression_transferred) / 1000
    experiments = float(experiments) / 1000
    experiments_transferred = float(experiments_transferred) / 1000
    database = float(database) / 1000
    database_transferred = float(database_transferred) / 1000
    textmining = float(textmining) / 1000
    textmining_transferred = float(textmining_transferred) / 1000
    initial_combined = int(initial_combined)


    ## compute prior away

    neighborhood_prior_corrected                 = compute_prior_away (neighborhood, prior)             
    neighborhood_transferred_prior_corrected     = compute_prior_away (neighborhood_transferred, prior) 
    fusion_prior_corrected                       = compute_prior_away (fusion, prior)             
    cooccurrence_prior_corrected                 = compute_prior_away (cooccurrence, prior)           
    coexpression_prior_corrected                 = compute_prior_away (coexpression, prior)            
    coexpression_transferred_prior_corrected     = compute_prior_away (coexpression_transferred, prior) 
    experiments_prior_corrected                  = compute_prior_away (experiments, prior)   
    experiments_transferred_prior_corrected      = compute_prior_away (experiments_transferred, prior) 
    database_prior_corrected                     = compute_prior_away (database, prior)      
    database_transferred_prior_corrected         = compute_prior_away (database_transferred, prior)
    textmining_prior_corrected                   = compute_prior_away (textmining, prior)            
    textmining_transferred_prior_corrected       = compute_prior_away (textmining_transferred, prior) 

    ## then, combine the direct and transferred scores for each category:

    neighborhood_both_prior_corrected = 1.0 - (1.0 - neighborhood_prior_corrected) * (1.0 - neighborhood_transferred_prior_corrected)
    coexpression_both_prior_corrected = 1.0 - (1.0 - coexpression_prior_corrected) * (1.0 - coexpression_transferred_prior_corrected)
    experiments_both_prior_corrected = 1.0 - (1.0 - experiments_prior_corrected) * (1.0 - experiments_transferred_prior_corrected)
    database_both_prior_corrected = 1.0 - (1.0 - database_prior_corrected) * (1.0 - database_transferred_prior_corrected)
    textmining_both_prior_corrected = 1.0 - (1.0 - textmining_prior_corrected) * (1.0 - textmining_transferred_prior_corrected)

    ## now, do the homology correction on cooccurrence and textmining:

    cooccurrence_prior_homology_corrected = cooccurrence_prior_corrected * (1.0 - homology)
    textmining_both_prior_homology_corrected = textmining_both_prior_corrected * (1.0 - homology)

    ## next, do the 1 - multiplication:

    combined_score_one_minus = (
    #    (1.0 - neighborhood_both_prior_corrected) *
    #    (1.0 - fusion_prior_corrected) *
    #    (1.0 - cooccurrence_prior_homology_corrected) *
        (1.0 - coexpression_both_prior_corrected) *
    #    (1.0 - experiments_both_prior_corrected) *
    #    (1.0 - database_both_prior_corrected) *
    #    (1.0 - textmining_both_prior_homology_corrected) 
	1) 

    ## and lastly, do the 1 - conversion again, and put back the prior *exactly once*

    combined_score = (1.0 - combined_score_one_minus)            ## 1- conversion
    combined_score *= (1.0 - prior)                              ## scale down
    combined_score += prior                                      ## and add prior.

    ## round

    combined_score = int(combined_score * 1000)
    #print(protein1, protein2, combined_score)
    list_to_print=(str(protein1), str(protein2), str(combined_score))
    file.write(','.join(list_to_print))
    file.write("\n")

 
file.close() 
 
