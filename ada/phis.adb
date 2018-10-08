with System.Storage_Elements;
with Custom_Types; use Custom_Types;

package body Phis is
   function Seed_To_Int ( Seed : Phi_Seed ) return Index is
   begin
      -- Seed is an array of bytes that is meant to be read in a big endian way
      return Index( (2**24)*Seed(1) + (2**16)*Seed(2) + (2**8)*Seed(3) + (2**0)*Seed(4) );
   end;

   function Phi ( Seed : Phi_Seed;
                  I : Index)
      -- this implementation of the Phi function comes from
      -- hhtps://www.cryptolux.org/images/0/0d/Argon2.pdf (page 7)
   return Index is
      R : Index := I-1;
      J : Index := Seed_To_Int ( Seed );
      X : Index := (J**2) / (2**32);
      Y : Index := (R*X) / (2**32);
   begin
      return R - Y;
   end;

end Phis;
