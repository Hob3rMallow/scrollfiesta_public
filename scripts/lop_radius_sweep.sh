#!/bin/bash
# Sweep MLS_PROJECT_RADIUS_VOX at fixed ITERS=50 to test the curvature-floor
# hypothesis: smaller R sees less curvature within the kernel so the local
# normal estimate is more accurate, projection direction is less biased,
# and final thickness should drop below the R=12 asymptote.
#
# The constant is on a single line "#define MLS_PROJECT_RADIUS_VOX  X.0f"
# (followed by a multi-line comment). We use sed to flip just the numeric
# value on that line.

cd "$(dirname "$0")/.."

OUT=output/lop_radius_sweep.txt
echo "=== LOP radius sweep (tangent projection, ITERS=50) ===" > $OUT
echo "Started: $(date)" >> $OUT
echo "" >> $OUT
echo "R     wall_s  pipe_s   comp001          comp002          comp003          comp004" >> $OUT

# Fix ITERS at 50 for this sweep
sed -i 's/^#define MLS_PROJECT_ITERS  *[0-9]\+$/#define MLS_PROJECT_ITERS       50/' src/common/pipeline_constants.h

for R in 3 5 7 9 11 13; do
    # Replace the radius value on the same line, preserve the comment
    sed -i "s/^#define MLS_PROJECT_RADIUS_VOX  [0-9.]*f /#define MLS_PROJECT_RADIUS_VOX  ${R}.0f /" src/common/pipeline_constants.h

    # Skip if R+1 > halo (would violate the cube/halo invariant)
    if [ "$R" -ge "14" ]; then
        echo "R=$R: skipping (would exceed halo=14)" | tee -a $OUT
        continue
    fi

    powershell.exe -NoProfile -Command "& 'C:\Program Files\Microsoft Visual Studio\2022\Professional\MSBuild\Current\Bin\MSBuild.exe' cube_mesh.vcxproj /p:Configuration=Release /p:Platform=x64 /m /v:minimal" > /tmp/build.log 2>&1 || {
        echo "R=$R: BUILD FAILED" | tee -a $OUT
        tail -10 /tmp/build.log >> $OUT
        continue
    }

    rm -rf output/lop_rad_$R
    mkdir -p output/lop_rad_$R
    TS=$(date +%s)
    ./build/Release/cube_mesh.exe PHerc0139-4x5x5/cubes_PRED/z04480_y03328_x02816.tif \
        output/lop_rad_$R/out.tif \
        --dump-obj output/lop_rad_$R/dump --no-qem --halo 14 > /tmp/run.log 2>&1
    TF=$(date +%s)
    ELAPSED=$((TF - TS))
    PIPE_TIME=$(grep -E "Timings" /tmp/run.log | tail -1 | awk '{for(i=1;i<=NF;i++) if($i~/total=/) print $i}' | sed 's/total=//; s/s$//')

    RESULTS=""
    for COMP in 001 002 003 004; do
        F="output/lop_rad_$R/dump/z04480_y03328_x02816/z04480_y03328_x02816_step0_post_lop_points/z04480_y03328_x02816_post_lop_points_comp$COMP.obj"
        if [ ! -f "$F" ]; then continue; fi
        awk '/^v / && $2>=4481.5 && $2<=4606.5 && $3>=3329.5 && $3<=3454.5 && $4>=2817.5 && $4<=2942.5' "$F" > /tmp/sw_$COMP.obj
        T=$(./scripts/lop_thickness_analyze.exe /tmp/sw_$COMP.obj 6.0 2>&1 | grep "Mean" | awk '{print $3}')
        BBOX=$(awk '/^v / && $2>=4481.5 && $2<=4606.5 && $3>=3329.5 && $3<=3454.5 && $4>=2817.5 && $4<=2942.5 {if(!s){zmin=$2;zmax=$2;ymin=$3;ymax=$3;xmin=$4;xmax=$4;s=1} if($2<zmin)zmin=$2; if($2>zmax)zmax=$2; if($3<ymin)ymin=$3; if($3>ymax)ymax=$3; if($4<xmin)xmin=$4; if($4>xmax)xmax=$4} END{printf "%.0fx%.0fx%.0f", zmax-zmin, ymax-ymin, xmax-xmin}' "$F")
        RESULTS="$RESULTS $T($BBOX)"
    done
    printf "R=%2d  %6d  %6s %s\n" $R $ELAPSED "$PIPE_TIME" "$RESULTS" | tee -a $OUT
done

echo "" >> $OUT
echo "Done: $(date)" >> $OUT
