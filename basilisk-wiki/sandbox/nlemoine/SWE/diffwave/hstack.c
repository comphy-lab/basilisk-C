int main()
{
  system ("wget http://basilisk.fr/sandbox/nlemoine/SWE/diffwave/fullsaintvenant/h.mp4 -O h_full.mp4 && "
          "wget http://basilisk.fr/sandbox/nlemoine/SWE/diffwave/diffwave/h.mp4 -O h_diffwave.mp4 && "
          "ffmpeg -i h_full.mp4 -i h_diffwave.mp4 -filter_complex hstack h_full_diffwave.mp4 > logffmpeg.txt &&"
          "rm h_full.mp4 h_diffwave.mp4");          
}

/**
![Video](hstack/logffmpeg.txt)
![Video](hstack/h_full_diffwave.mp4)(width=80% )
*/