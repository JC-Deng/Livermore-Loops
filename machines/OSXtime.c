 double (*myNVtime)();
         SV **svp = hv_fetch(PL_modglobal, "Time::NVtime", 12, 0);
         if (!svp)         croak("Time::HiRes is required");
         if (!SvIOK(*svp)) croak("Time::NVtime isn't a function pointer");
         myNVtime = INT2PTR(double(*)(), SvIV(*svp));
         printf("The current time is: %f\n", (*myNVtime)());

