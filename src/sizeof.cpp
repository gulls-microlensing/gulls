/* A program to print the size of the different    */
/* kinds of C variables on your machine.           */
#include <stdio.h>
#include <limits.h>

int main(void)
{   
    printf ("Signed   : Size %20s %22s\n", "Min", "Max") ;
    printf ("char     :  %d %22d %22d\n", 
            (int) sizeof (char), CHAR_MIN,CHAR_MAX); 
    printf ("short    :  %d %22d %22d\n", 
            (int) sizeof (short), SHRT_MIN,SHRT_MAX); 
    printf ("int      :  %d %22d %22d\n", 
            (int) sizeof (int), INT_MIN,INT_MAX); 
    printf ("long     :  %d %22ld %22ld\n", 
            (int) sizeof (long), LONG_MIN,LONG_MAX); 
    printf ("\n") ;
    
    printf ("Unsigned : Size %20s %22s\n", "Min", "Max") ;
    printf ("char     :  %d %22d %22u\n", 
            (int) sizeof (unsigned char),0,UCHAR_MAX); 
    printf ("short    :  %d %22d %22u\n", 
            (int) sizeof (unsigned short),0,USHRT_MAX); 
    printf ("int      :  %d %22d %22u\n", 
            (int) sizeof (unsigned int),0,UINT_MAX); 
    printf ("long     :  %d %22d %22lu\n", 
            (int) sizeof (unsigned long),0,ULONG_MAX); 
    printf ("\n") ;

    printf ("single prec. float : %d\n", (int) sizeof (float)); 
    printf ("double prec. float : %d\n", (int) sizeof (double)); 

    return 0 ;
}
