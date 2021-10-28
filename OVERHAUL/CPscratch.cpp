if (pj.e_sub>efcheck)	// impose freeze-out check for e, not s
    {
      pj.Freeze=0;
    }
    else
    {
      pj.Freeze=4;
      --kk;
      ++numpart;
    }