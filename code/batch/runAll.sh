rm /home/livingstonb/GitHub/Continuous_Time_HA/output/two_asset/*
rm /home/livingstonb/GitHub/Continuous_Time_HA/output/con_effort/*

sbatch /home/livingstonb/GitHub/Continuous_Time_HA/code/two_asset/table_tests.sbatch
sbatch /home/livingstonb/GitHub/Continuous_Time_HA/code/con_effort/server.sbatch
