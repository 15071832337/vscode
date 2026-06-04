#!bin/bash

export OPTIONS="-b --configuration json://${PWD}/configuration_sourav_unweight_new.json --resources-monitoring 2 --aod-memory-rate-limit 1000000000 --shm-segment-size 6000000000"
o2-analysis-lf-lambdaspincorrderived ${OPTIONS}  > --aod-file @input_data.txt >stdout.log 2>&1 --aod-writer-json OutputDirector.json



















#--resources-monitoring 2 --aod-memory-rate-limit 1000000000 --shm-segment-size 6000000000
#|o2-analysis-pid-tof-merge ${OPTIONS} \
#| o2-analysis-pid-tof ${OPTIONS} \
#| o2-analysis-pid-tof-full ${OPTIONS} \
#| o2-analysis-pid-tof-base ${OPTIONS} \
#| o2-analysis-lf-v0qaanalysis ${OPTIONS} \
#| o2-analysis-lf-v0postprocessing ${OPTIONS} \
#| o2-analysis-lf-strangenesstofpid ${OPTIONS} \
#| o2-analysis-lf-strangepidqa ${OPTIONS} \
#-aod-file @input1.txt >stdout.log 2>&1 --aod-writer-keep AOD/MYTABLE/0,AOD/MYTABLEJet/0,AOD/LeadingJet/0
#o2-analysis-propagationservice ${OPTIONS} \
#| o2-analysis-multiplicity-table ${OPTIONS} \
#| o2-analysis-trackselection ${OPTIONS} \
#| o2-analysis-event-selection-service ${OPTIONS} \
#| o2-analysis-centrality-table ${OPTIONS} \
#| o2-analysis-tracks-extra-v002-converter ${OPTIONS} \
#| o2-analysis-lf-lambdajetpolarizationbuilder ${OPTIONS} \
#| o2-analysis-lf-strangenesstofpid ${OPTIONS} \
#| o2-analysis-lf-strangepidqa ${OPTIONS} \
#| o2-analysis-lf-strangederivedbuilder ${OPTIONS} \
#| o2-analysis-pid-tpc-service ${OPTIONS} \
#| o2-analysis-pid-tof-merge ${OPTIONS} \
#| o2-analysis-ft0-corrected-table ${OPTIONS} \
#| o2-analysis-je-emcal-correction-task  ${OPTIONS} \
#| o2-analysis-je-jet-deriveddata-producer ${OPTIONS} \
#| o2-analysis-je-jet-finder-data-charged ${OPTIONS} \
#| o2-analysis-je-jet-finder-mcp-charged ${OPTIONS} \
#| o2-analysis-je-jet-finder-mcd-charged ${OPTIONS} \
#| o2-analysis-mm-track-propagation ${OPTIONS} \
#| o2-analysis-track-to-collision-associator ${OPTIONS} \
#| o2-analysis-fwdtrack-to-collision-associator ${OPTIONS} \
#| o2-analysis-lf-lambdajetpolarization ${OPTIONS}> 








#| o2-analysis-pid-tpc-base ${OPTIONS} \
#| o2-analysis-pid-tof-base ${OPTIONS} \
#| o2-analysis-pid-tof ${OPTIONS} \
#| o2-analysis-pid-tof-full ${OPTIONS} \
#| o2-analysis-je-emcal-correction-task  ${OPTIONS} \
#| o2-analysis-pid-tpc ${OPTIONS} \
#| o2-analysis-lf-lambdakzerobuilder ${OPTIONS} \
#o2-analysis-timestamp ${OPTIONS} \
#| o2-analysis-event-selection ${OPTIONS} \
#o2-analysis-track-propagation ${OPTIONS} \
#o2-analysistutorial-mm-my-example-task --aod-file AO2D1.root
#/o2-analysis-timestamp ${OPTIONS} \
#| o2-analysis-tracks-extra-v002-converter ${OPTIONS} \
#| o2-analysis-track-propagation ${OPTIONS} \
#| o2-analysis-je-jet-matching-mc ${OPTIONS} \
#| --aod-writer-keep AOD/MYTABLE/0
#| o2-analysis-multiplicity-table ${OPTIONS} \
#| o2-analysis-lf-lambdajetpolarization
#| o2-analysistutorial-lf-myanalysis
#| o2-analysis-je-emcal-correction-task  ${OPTIONS} \
#| o2-analysis-je-jet-deriveddata-producer ${OPTIONS} \
#| o2-analysis-je-jet-finder-data-charged ${OPTIONS} \
#| o2-analysis-je-jet-finder-mcp-charged ${OPTIONS} \
#| o2-analysis-je-jet-finder-mcd-charged ${OPTIONS} \
#| o2-analysis-je-jet-matching-mc ${OPTIONS} \

#| o2-analysis-je-jet-finder-data-full ${OPTIONS} \
#| o2-analysis-je-jet-tutorial ${OPTIONS} \
#| o2-analysis-je-jet-tutorial ${OPTIONS} \
#| o2-analysis-je-phi-in-jets ${OPTIONS} \
#| o2-analysis-je-jet-finder-charged-qa ${OPTIONS}\
#| o2-analysis-je-phi-in-jets ${OPTIONS} \

#| o2-analysis-je-phi-in-jets ${OPTIONS} \
#| o2-analysis-pid-tpc-full ${OPTIONS} \
#| o2-analysis-v0converter ${OPTIONS} \
#| o2-analysis-bc-converter ${OPTIONS} \
#| o2-analysis-tracks-extra-converter ${OPTIONS} \


#| o2-analysis-v0converter ${OPTIONS} \
#| o2-analysis-centrality-table ${OPTIONS}
#| o2-analysis-je-jet-deriveddata-producer ${OPTIONS} \
#| o2-analysis-lf-derivedlambdakzeroanalysis ${OPTIONS} \
#| o2-analysis-lf-lambdakzeropid ${OPTIONS} \
#| o2-analysis-pid-tof-base ${OPTIONS} \
#| o2-analysis-pid-tpc-full ${OPTIONS} \
#| o2-analysis-lf-strangederivedbuilder ${OPTIONS} \
#| o2-analysis-lf-cascadebuilder ${OPTIONS} \
#| o2-analysis-lf-cascademcbuilder ${OPTIONS} \
#| o2-analysis-lf-lambdakzeromcbuilder ${OPTIONS} \
#| o2-analysistutorial-lf-strangeness-mytask-1 ${OPTIONS} \
#| o2-analysis-lf-cascademcbuilder ${OPTIONS} \
#| o2-analysis-lf-lambdakzeropid ${OPTIONS} \
#| o2-analysis-pid-tof-base ${OPTIONS} \
#| o2-analysis-pid-tpc-full ${OPTIONS} \
#| o2-analysis-lf-strangederivedbuilder ${OPTIONS} \
#| o2-analysis-lf-cascadebuilder ${OPTIONS} \
#| o2-analysis-lf-cascademcbuilder ${OPTIONS} \
#| o2-analysis-lf-lambdakzeromcbuilder ${OPTIONS} \
#| o2-analysis-lf-derivedlambdakzeroanalysis ${OPTIONS}

#| o2-analysis-tracks-extra-converter ${OPTIONS} \
#| o2-analysis-lf-lambdakzeromcbuilder ${OPTIONS} \
#| o2-analysis-je-jet-finder-data-charged ${OPTIONS} \
#| o2-analysis-je-emcal-correction-task ${OPTIONS} \
#| o2-analysis-je-jet-finder-data-charged ${OPTIONS} \
#| o2-analysis-trackselection ${OPTIONS} \
#| o2-analysis-centrality-table ${OPTIONS} \
#| o2-analysis-je-jet-deriveddata-producer ${OPTIONS}>stdout.log 2>&1
#| o2-analysis-zdc-converter ${OPTIONS}\
#| o2-analysis-lf-lambdakzeropid ${OPTIONS} \
#| o2-analysis-pid-tpc-full ${OPTIONS} \
#| o2-analysis-pid-tof-base ${OPTIONS} \
#| o2-analysis-lf-strangederivedbuilder ${OPTIONS} \
#| o2-analysis-lf-lambdakzeromcbuilder ${OPTIONS} \
#| o2-analysis-lf-cascademcbuilder ${OPTIONS} \
#| o2-analysis-zdc-converter ${OPTIONS}\
#| o2-analysis-centrality-table ${OPTIONS}> stdout.log 2>&1
#|o2-analysis-pid-tpc-full ${OPTIONS} \
#| o2-analysis-pid-tof-base ${OPTIONS} \
#| o2-analysis-centrality-table ${OPTIONS}> stdout.log 2>&1
#| o2-analysis-lf-lambdakzeropid ${OPTIONS} \
#| o2-analysis-pid-tpc-full ${OPTIONS} \
#| o2-analysis-pid-tof-base ${OPTIONS} \
#| o2-analysis-lf-strangederivedbuilder ${OPTIONS} \
#| o2-analysis-lf-lambdakzeromcbuilder ${OPTIONS} \
#| o2-analysis-lf-cascademcbuilder ${OPTIONS} > stdout.log 2>&1
#| o2-analysis-centrality-table ${OPTIONS} > stdout.log 2>&1
#| o2-analysis-centrality-table ${OPTIONS}> stdout.log 2>&1

#| o2-analysis-lf-strangederivedbuilder ${OPTIONS}\
#| o2-analysis-pid-tpc-full ${OPTIONS} \
#| o2-analysis-centrality-table ${OPTIONS}>stdout.log 2>&1

#| o2-analysis-lf-strangederivedbuilder ${OPTIONS} \
#| o2-analysis-pid-tpc-full ${OPTIONS} \
#| o2-analysis-pid-tof-base ${OPTIONS} \
#| o2-analysis-lf-lambdakzeromcbuilder ${OPTIONS} \
#| o2-analysis-lf-cascademcbuilder ${OPTIONS}\
#| o2-analysis-centrality-table ${OPTIONS}>stdout.log 2>&1


#| o2-analysis-lf-cascadebuilder ${OPTIONS} \
#| o2-analysis-centrality-table ${OPTIONS}> stdout.log 2>&1

#| o2-analysis-pid-tpc ${OPTIONS}

