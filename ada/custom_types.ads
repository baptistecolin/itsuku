with Interfaces;

package Custom_Types is
   type Index is new Natural;
   type Byte is new Interfaces.Unsigned_8;
   type Byte_Array is array (Natural range <>) of Byte;
end Custom_Types;
