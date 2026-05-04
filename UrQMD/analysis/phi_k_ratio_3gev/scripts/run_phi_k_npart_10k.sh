#!/usr/bin/env sh
set -eu

run_dir="${1:-/home/starsdong.guest/runs/urqmd-4.0_3gev_1k_f16_phi}"
cd "$run_dir"

rm -f test.f13 test.f14 test.f15 test.f16 test.f19 test.f20 run_npart_batch_*.log inputfile phi_k_npart_events.tsv
sync
sudo fstrim -av >/dev/null 2>&1 || true

printf "batch\tseed\tevent\tb_fm\tnpart\tphi\tKminus\tKplus\tK0\tantiK0\n" > phi_k_npart_events.tsv

for batch in 1 2 3 4 5 6 7 8 9 10; do
  seed=$((246800 + batch * 9973))
  rm -f test.f16 inputfile "run_npart_batch_${batch}.log"
  {
    printf "pro 197 79\n"
    printf "tar 197 79\n\n"
    printf "nev 1000\n"
    printf "imp -15.0\n"
    printf "ecm 3.0\n"
    printf "tim 200 200\n"
    printf "eos 0\n"
    printf "rsd %d\n\n" "$seed"
    printf "f13\nf14\nf15\nf19\nf20\n\n"
    printf "xxx\n"
  } > inputfile

  printf "starting batch %d seed %d\n" "$batch" "$seed"
  ./runqmd.bash > "run_npart_batch_${batch}.log" 2>&1

  awk -v batch="$batch" -v seed="$seed" '
    function emit() {
      if (have_event) {
        printf "%d\t%d\t%d\t%.6f\t%.1f\t%d\t%d\t%d\t%d\t%d\n",
          batch, seed, event, b, npart, phi, km, kp, k0, ak0
      }
      phi = km = kp = k0 = ak0 = 0
    }
    /impact_parameter_real/ {
      emit()
      b = $2 + 0
      npart = 0
      have_event = 0
      next
    }
    /Participants_Glauber/ {
      npart = $3 + 0
      next
    }
    /event#/ {
      event = $2 + 0
      have_event = 1
      next
    }
    NF >= 15 && $10 ~ /^-?[0-9]+$/ && $11 ~ /^-?[0-9]+$/ {
      if ($10 == 109 && $11 == 0) phi++
      if ($10 == -106 && $11 == -1) km++
      if ($10 == 106 && $11 == 1) kp++
      if ($10 == 106 && $11 == -1) k0++
      if ($10 == -106 && $11 == 1) ak0++
    }
    END {
      emit()
    }
  ' test.f16 >> phi_k_npart_events.tsv

  awk -F '\t' -v batch="$batch" '
    NR > 1 && $1 == batch {
      ev++
      phi += $6
      km += $7
      np += $5
    }
    END {
      printf "batch %d: events=%d <Npart>=%.2f phi=%d Kminus=%d ratio=%.6f\n",
        batch, ev, np / ev, phi, km, (km ? phi / km : 0)
    }
  ' phi_k_npart_events.tsv

  rm -f test.f16 "run_npart_batch_${batch}.log" inputfile
  sync
  sudo fstrim -av >/dev/null 2>&1 || true
done

printf "done\n"
