package de.labathome;

import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;
import java.lang.annotation.ElementType;

/**
 * Annotation used to identify namelist variables for parsing definition classes.
 * In Fortran, indices may start at any(?) arbitrary integer number, i.e. also at negative numbers.
 * The starting indices can be marked by the dim0min(for 1d and 2d arrays) and dim1min (for 2d arrays).
 * These should also be used if an array is indexed beginning with 0.
 * @author Jonathan Schilling (jonathan.schilling@ipp.mpg.de)
 * @version 0.9 2018-07-25 intial implementation
 * @version 1.0 2018-09-11 1d and 2d arrays with fancy specifiers
 */
@Retention(RetentionPolicy.RUNTIME)
@Target(ElementType.FIELD)
public @interface namelist_variable {

    String name() default "[variable_name]";
    int dim0min() default 1; // starting index for dimension 0, e.g. row
    int dim1min() default 1; // starting index for dimension 1, e.g. col
}
