for i in {0..16}
  do
		./absorptionPixelSum Cu 1 1 1 1 $i 1
		./absorptionPixelSum Fe 1 1 1 1 $i 1
		./absorptionPixelSum Al 1 1 1 1 $i 1
		./absorptionPixelSum Sand 1 1 1 1 $i 1
		./absorptionPixelSum H2O 1 1 1 1 $i 1
		./absorptionPixelSum Paraffin 1 1 1 1 $i 1
 done

