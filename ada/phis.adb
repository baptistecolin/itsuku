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
   return Index is
      R : Index := I-1;
   begin
      return I-1; -- TODO : write actual body
   end;

end Phis;
