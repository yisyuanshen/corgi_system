# /bin/bash
trap "kill -- -$$" INT

export CORE_LOCAL_IP=127.0.0.1
export CORE_MASTER_ADDR=127.0.0.1:10010

NodeCore &
sleep .1
echo Enter to continue...
read
pushd build
./FSM_node_webot &
sleep .1
./FSM_listen_node_webot &
sleep 1
# # ./state_estimation_webot &
# # sleep .1
popd
pushd ~/corgi_ws/corgi_system/js/control_panel
node quadruped.js &
popd

# pushd build
# ./Kuramoto_test_webot &
# sleep .1
# ./state_kf_webot &
# sleep .1
# popd

# pushd build
# ./my_cpg_test &
# sleep .1
# ./state_kf_webot &
# sleep .1
# popd

wait
