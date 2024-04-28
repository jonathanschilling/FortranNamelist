package de.labathome;

import java.lang.annotation.Annotation;
import java.lang.reflect.Field;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Locale;
import java.util.Objects;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Class to parse Fortran namelists into Java classes.
 * This uses the information given in the class definition to constrain the parser.
 * @author Jonathan Schilling (jonathan.schilling@mail.de)
 * @version 0.9.0 2018-07-25 intial implementation
 * @version 1.0.1 2018-09-11 1d and 2d arrays with fancy specifiers
 */
public class FortranNamelist {

    private String namelist;
    private String groupName;
    private Object parseInto;

    /**
     * set to true to be enable verbal diarrhea with lots of debug output
     */
    public boolean _debug = false;

    /**
     * Define parser for group {@code groupName} from the namelist given as text in {@code namelist}
     * according to the Object definition given as {@code parseInto}.
     * @param _namelist String containing the namelist
     * @param _groupName group name of the group to be parsed into the object {@code parseInto}
     * @param _parseInto Object definition for parsing target. Use annotation {@code namelist_variable(name="lstell_sym")}
     *                  to specify which variables to parse. {@code name} is optional (default=take variable name as identifier
     *                  in namelist) and can be used to map from Fortran names (e.g. lstell_sym) to your own defitions (e.g. stellaratorSymmetric).
     */
    public FortranNamelist(String _namelist, String _groupName, Object _parseInto) {
        this.namelist = _namelist;
        this.groupName = _groupName.toLowerCase();
        this.parseInto = _parseInto;
    }

    /**
     * Call the parser and return the parsed Object.
     * @return Object of same type as {@code parseInto} with member variables filled out as stated in namelist.
     */
    public Object getParsed() {
        parseIt();
        return parseInto;
    }

    /**
     * Actually parse the group {@code groupName} from the namelist given in {@code namelist} into an Object of type {@code parseInto}.
     */
    private void parseIt() {
        // only work if all input is given
        if (namelist == null || namelist.equals("") || groupName == null || groupName.equals("") || parseInto == null) {
            throw new RuntimeException("FortranNamelist needs namelist and groupName and parseInto to be set!");
        }

        // map of variable names to look for in the namelist
        // key is variable name in the namelist and value is field name in parseInto
        HashMap<String, String> names = new HashMap<>();

        // map of variable types to look for in the namelist
        // key is variable name in the namelist and value is name of type in parseInto
        HashMap<String, String> types = new HashMap<>();

        HashMap<String, Integer> dim0min = new HashMap<>();
        HashMap<String, Integer> dim1min = new HashMap<>();

        // first step: create definition of what to look for in the namelist
        for (Field field : parseInto.getClass().getDeclaredFields()) {
            // request access to all member fields
            field.setAccessible(true);

            try {
                // loop over annotations of field
                for (Annotation annotation : field.getAnnotations()) {
                    // check for namelist_variable annotation
                    if (annotation instanceof namelist_variable) {
                        namelist_variable nmlAnnotation = (namelist_variable) annotation;

                        // get name and check wether variable name is to be used
                        String name = nmlAnnotation.name().toLowerCase();
                        String defaultName = (String) namelist_variable.class.getDeclaredMethod("name").getDefaultValue();
                        if (name.equals(defaultName.toLowerCase())) {
                            // System.out.println("Field "+field.getName()+" is called \"" + field.getName().toLowerCase() + "\" in namelist.");
                            name=field.getName().toLowerCase();
                        }
                        //System.out.print(name + " => " + field.getName()+" is ");
                        names.put(name, field.getName());

                        // Starting indices for >0d arrays, since Fortran allows for negative indices...
                        int dim0min_val = nmlAnnotation.dim0min();
                        int dim1min_val = nmlAnnotation.dim1min();
                        dim0min.put(name, dim0min_val);
                        dim1min.put(name, dim1min_val);

                        // select which target type to try to fill from the namelist contents
                        if        (field.getType().equals(int.class)) {
                            types.put(name, "0d:int");
                        } else if (field.getType().equals(int[].class)) {
                            types.put(name, "1d:int");
                        } else if (field.getType().equals(int[][].class)) {
                            types.put(name, "2d:int");
                        } else if (field.getType().equals(double.class)) {
                            types.put(name, "0d:double");
                        } else if (field.getType().equals(double[].class)) {
                            types.put(name, "1d:double");
                        } else if (field.getType().equals(double[][].class)) {
                            types.put(name, "2d:double");
                        } else if (field.getType().equals(String.class)) {
                            types.put(name, "0d:String");
                        } else if (field.getType().equals(String[].class)) {
                            types.put(name, "1d:String");
                        } else if (field.getType().equals(String[][].class)) {
                            types.put(name, "2d:String");
                        } else if (field.getType().equals(boolean.class)) {
                            types.put(name, "0d:boolean");
                        } else if (field.getType().equals(boolean[].class)) {
                            types.put(name, "1d:boolean");
                        } else if (field.getType().equals(boolean[][].class)) {
                            types.put(name, "2d:boolean");
                        } else {
                            System.out.println("Type of \""+ name + "\" was not recognized: " + field.getAnnotatedType());
                        }
                    }
                }
            } catch (Exception e) {
            	throw new RuntimeException(e);
            }
        }

        if (_debug) {
            System.out.println("Variable names to look for:");
            System.out.println(names.keySet());

            System.out.println("Types to parse them into:");
            System.out.println(types.values());
        }

        // remove single-line comments like ! this is a comment or
        // ! this is a comment without stuff in front of it
        while (namelist.indexOf("!") > -1) {
            int commentStart = namelist.indexOf("!");
            int lineEnd = namelist.indexOf("\n", commentStart);
            if (lineEnd == -1) {
            	// no newline is at the end of the file
            	lineEnd = namelist.length() - 1;
            }

            if (_debug) {
                System.out.println("INFO: found !-comment from " + commentStart + " to " + lineEnd+":");
                System.out.println(namelist.substring(commentStart, lineEnd+1));
            }

            // remove comment from String inputFile
            String commentRemoved = namelist.substring(0, commentStart) + namelist.substring(lineEnd+1);
            namelist = commentRemoved;
        }

        // remove single-line comments like
        // c this is a comment without stuff in front of it
        while (namelist.indexOf("\nc ") > -1 || namelist.indexOf("\nc\r\n") > -1 || namelist.indexOf("\nc\n") > -1) {
            int commentStart = namelist.indexOf("\nc\n");
            int lineEnd = commentStart+1;
            if (commentStart < 0) {
                commentStart = namelist.indexOf("\nc\r\n");
                lineEnd = commentStart+2;
            }
            if (commentStart < 0) {
                commentStart = namelist.indexOf("\nc ")+1;
                lineEnd = namelist.indexOf("\n", commentStart);
                if (lineEnd == -1) {
                	// no newline is at the end of the file
                	lineEnd = namelist.length() - 1;
                }
            }

            if (_debug) {
                System.out.println("INFO: found c-comment from " + commentStart + " to " + lineEnd+":");
                System.out.println(namelist.substring(commentStart, lineEnd+1));
            }

            // remove comment from String inputFile
            String commentRemoved = namelist.substring(0, commentStart) + namelist.substring(lineEnd+1);
            namelist = commentRemoved;
        }

        // replace " by ' for unified parsing
        namelist = namelist.replace("\"", "'");

        int namelistStart = 0, nameEnd = -1;
        boolean foundNamelist = false;
        // loop over all namelists in the file and interpret only the groupName one
        while (!foundNamelist && namelistStart < namelist.length()) {
	        // look for starting tags of namelists like "&input \n"
	        namelistStart = namelist.indexOf("&", namelistStart)+1;
	        nameEnd = namelist.indexOf("\n", namelistStart);
	        // check wether accidentially a variable name was included (variables may be already defined in the first line,
	        // right after the namelist starting tag
	        if (namelist.substring(namelistStart, nameEnd).contains(" ")) {
	            nameEnd = namelist.indexOf(" ", namelistStart);
	        }

	        // get group name
	        String namelistName = namelist.substring(namelistStart, nameEnd).trim().toLowerCase();

	        if (_debug) {
	        	System.out.println("found namelist: " + namelistName);
	        }

	        if (namelistName.equals(groupName)) {
	        	foundNamelist = true;
	        } else {
	        	// skip the start of the current namelist to start searching for the next namelist in the input
	        	namelistStart++;
	        }
        }

        if (foundNamelist) {
            // find end of namelist: "/"
            int namelistEnd = nameEnd;
            boolean insideString = false;
            while (namelistEnd < namelist.length()) {
                String currentChar = namelist.substring(namelistEnd, namelistEnd+1);
                if (currentChar.equals("\'")) {
                    insideString = !insideString;
                }
                if (!insideString && currentChar.equals("/")) {
                    break;
                } else {
                    namelistEnd++;
                }
            }

            String contents = namelist.substring(nameEnd, namelistEnd);

            // God praise https://regex101.com/ !

            // look for a floating point number with exponential, p.ex. -1.2e-3 or 123.45E8 1.d-10
            Pattern floatingPointNumberPattern = Pattern.compile("^[-+]?[0-9]*(\\.)?[0-9]+([dDeE][-+]?[0-9]+)?");

            // in a String "1234.234 myvar", look for the varname "myvar"
            // and also array index varnames, p.ex. "myarr(-2, 4:6)"
//          Pattern nextVarNamePattern = Pattern.compile("([a-z_]+[0-9_]*[a-z_]*)\\s*(\\((\\s*[-+]?\\s*[0-9]+(\\s*:(\\s*[-+]?\\s*[0-9]+))*)?\\s*(,?(\\s*[-+]?\\s*[0-9]+(\\s*:(\\s*[-+]?\\s*[0-9]+))*))\\s*\\))?$");
            Pattern nextVarNamePattern = Pattern.compile("([a-z_]+[0-9_]*[a-z_]*)\\s*(\\((\\s*[-+]?\\s*[0-9]*(\\s*:(\\s*[-+]?\\s*[0-9]*))*)?\\s*(,(\\s*[-+]?\\s*[0-9]*(\\s*:(\\s*[-+]?\\s*[0-9]*))*)?)*\\s*\\))?$");

            int i=0, last_i=0;
            int var_counter=0;
            String name="", val="";

            // loop over file contents
            // Inside this loop, the content is split at "=" signs.
            // The first part is interpreted as a variable name.
            // Each next part is then split into a value of the variable from the iteration before and the next variable name.
            // The last part is interpreted as the value of the last variable.
            while (i<contents.length()) {

                last_i=i;
                // loop to next "="
                while (i<contents.length() && contents.charAt(i) != '=') i++;

                // at the beginning, look at the name of the first variable
                if (var_counter==0) {
                    name=contents.substring(last_i, i).trim().toLowerCase();
                } else {
                    // next part to split into some variable content (belonging to the variable identified by the current value of "name")
                    // and the next variable name
                    String part = contents.substring(last_i, i).trim();

                    // look for String definition like 'a String'
                    if (part.contains("'")) {
                        int begin = part.indexOf("'");
                        int end   = part.indexOf("'", begin+1);

                        String stringPart = part.substring(begin+1, end);

                        // now check wether there was something between the string start and the = sign
                        if (part.substring(0, begin).trim().length() > 0) {
                            System.out.println("Something is wrong, why is there something between a \"=\" and a string start in the input file ?" + part.substring(0, begin) );
                        } else {
                            val = stringPart.trim();

                            if (_debug) System.out.println(name + " => \"" + val + "\"");

                            // only try to interpret variables that are in the parseInto class
                            if (names.containsKey(name)) {
                                String fieldName = names.get(name);
                                String type = types.get(name);
                                if (_debug) System.out.println("interpret " + name + " as " + type);

                                // check for correct type (so whether we are actually looking for a string)
                                if (type.equals("0d:String")) {

                                    // set field to value from namelist
                                    try {
                                        Field field = parseInto.getClass().getDeclaredField(fieldName);
                                        field.setAccessible(true);
                                        field.set(parseInto, val);
                                    }
                                    catch (Exception e) {
                                        e.printStackTrace();
                                    }
                                }
                            }

                            // next variable name
                            if (end<part.length()) {
                                name = part.substring(end+1).trim().toLowerCase();
                                // some namelists have their entries separated by ",", so remove them in this ste
                                if (name.startsWith(",")) {
                                    name = name.substring(1).trim();
                                }
                                if (_debug) System.out.println("next name: " + name);
                            }
                        }
                    }

                    // look for logical definition like ".true." or ".false."
                    else if (part.toLowerCase().startsWith(".t") || part.toLowerCase().startsWith(".f")) {
                        int end = part.indexOf(".", 1);

                        String logicalPart = part.substring(0, end+1);

                        val = logicalPart;

                        if (_debug) System.out.println(name + " => \"" + val + "\"");

                        // only try to interpret variables that are in the parseInto class
                        if (names.containsKey(name)) {
                            String fieldName = names.get(name);
                            String type = types.get(name);
                            if (_debug) System.out.println("interpret " + name + " as " + type);

                            // check for correct type (so whether we are actually looking for a string)
                            if (type.equals("0d:boolean")) {

                                boolean value = (val.toLowerCase().startsWith(".t") ? true : false);

                                // set field to value from namelist
                                try {
                                    Field field = parseInto.getClass().getDeclaredField(fieldName);
                                    field.setAccessible(true);
                                    field.set(parseInto, value);
                                }
                                catch (Exception e) {
                                    e.printStackTrace();
                                }
                            }
                        }

                        // next variable name
                        if (end<part.length()) {
                            name = part.substring(end+1).trim().toLowerCase();
                            // some namelists have their entries separated by ",", so remove them in this ste
                            if (name.startsWith(",")) {
                                name = name.substring(1).trim();
                            }
                            if (_debug) System.out.println("next name: " + name);
                        }
                    }

                    // look for single-character logicals like "t" or "f"
                    else if (part.toLowerCase().startsWith("t") || part.toLowerCase().startsWith("f")) {
                        val = part.toLowerCase().substring(0, 1);

                        if (_debug) System.out.println(name + " => \"" + val + "\"");

                        // only try to interpret variables that are in the parseInto class
                        if (names.containsKey(name)) {
                            String fieldName = names.get(name);
                            String type = types.get(name);
                            if (_debug) System.out.println("interpret " + name + " as " + type);

                            // check for correct type (so whether we are actually looking for a string)
                            if (type.equals("0d:boolean")) {
                                boolean value = (val.toLowerCase().startsWith("t") ? true : false);

                                // set field to value from namelist
                                try {
                                    Field field = parseInto.getClass().getDeclaredField(fieldName);
                                    field.setAccessible(true);
                                    field.set(parseInto, value);
                                }
                                catch (Exception e) {
                                    e.printStackTrace();
                                }
                            }
                        }

                        // next variable name
                        name = part.substring(1).trim().toLowerCase();
                        // some namelists have their entries separated by ",", so remove them in this ste
                        if (name.startsWith(",")) {
                            name = name.substring(1).trim();
                        }
                        if (_debug) System.out.println("next name: " + name);

                    }

                    // look for (floating point) numbers
                    else if (floatingPointNumberPattern.matcher(part.toLowerCase()).find()) {

                        part=part.replace("\n", " ");
                        part=part.replace("\r", " ");

                        if (_debug) System.out.println("found number(s) as value: " + part.toLowerCase());

                        // it may be better to look for non-number character (arrays...) first...

                        int end=part.length();

                        // look using regex for next var name
                        Matcher m = nextVarNamePattern.matcher(part.toLowerCase());

                        if (m.find()) {
                            end = m.start()-1;
                            if (_debug) System.out.println("found next var name: " + part.substring(end));
                        }

                        if (end<0) {
                            end = part.length();
                        }
                        val = part.toLowerCase().substring(0, end).trim();

                        // 1.d-10 --> 1.e-10
                        val = val.replace("d", "e");
                        val = val.replace("D", "e");

                        // remove possible "," at end of strings, if someone might have felt a subtle need for separators
                        if (val.endsWith(",")) {
                            val = val.substring(0, val.length()-1).trim();
                        }

                        // remove possible ";" at end of strings, if someone might have felt a subtle need for separators
                        if (val.endsWith(";")) {
                            val = val.substring(0, val.length()-1).trim();
                        }

                        if (_debug) System.out.println(name + " => \"" + val + "\"");

                        // look for multi-index notation like 1.0 1.5 2*2.0 1.0 for later expansion
                        //                                            ^^^^^
                        // into normal array format 1.0 1.5 2.0 2.0 1.0
                        //                                  ^^^^^^^
                        boolean multiIndex = val.contains("*");

                        // look ahead for array element specifiers like myarr(6,5)=3.4 or yourarr(1:4)=10.0
                        boolean arrayIndexSpecifierFound = name.contains("(") && name.substring(name.indexOf("(")+1).contains(")");
                        String indexString = "";
                        if (arrayIndexSpecifierFound) {
                            indexString = name.substring(name.indexOf("(")+1, name.indexOf(")"));
                            name = name.substring(0, name.indexOf("(")).trim();
                            if (_debug) System.out.println(" ==> index in array \""+indexString+"\"");
                        }

                        // only try to interpret variables that are in the parseInto class
                        if (names.containsKey(name)) {
                            String fieldName = names.get(name);
                            String type = types.get(name);
                            if (_debug) System.out.println("interpret " + name + " as " + type);

                            try {
                                // check for correct type and parse the value(s) accordingly
                                if (type.equals("0d:int")) {
                                    int value = Integer.valueOf(val);

                                    // set field to value from namelist
                                    Field field = parseInto.getClass().getDeclaredField(fieldName);
                                    field.setAccessible(true);
                                    field.set(parseInto, value);

                                } else if (type.equals("0d:double")) {
                                    double value = Double.valueOf(val);

                                    // set field to value from namelist
                                    Field field = parseInto.getClass().getDeclaredField(fieldName);
                                    field.setAccessible(true);
                                    field.set(parseInto, value);

                                } else if (type.equals("1d:int")) {
                                    int dim0min_val = dim0min.get(name);

                                    // get field to value from namelist
                                    Field field = parseInto.getClass().getDeclaredField(fieldName);
                                    field.setAccessible(true);
                                    int[] value = (int[])field.get(parseInto);

                                    // two possibilities:
                                    if (arrayIndexSpecifierFound) {
                                        // a1) indexed value: myarr(  3) = 5.64e4  => insert at given position
                                        // a2)   range index: myarr(3:5) = 0.0     => insert into all given positions
                                        if (indexString.contains(",")) {
                                            if (_debug) System.out.println("error: try to put value into 1d "+name+" at location indicated by more than one index: " + indexString);
                                        } else {
                                            // check for range specifiers in indices
                                            if (indexString.contains(":")) {
                                            	String[] inValues = val.split("[\\s,]+");
                                            	int idxStart;
                                            	int idxEnd;
                                            	if (indexString.length() == 1) {
                                            		// case a3) wildcard: "myarr(:) = 1.0, 2.0, 3.0"
                                            		idxStart = dim0min_val;
                                            		idxEnd = inValues.length-1 + dim0min_val;
                                            	} else {
                                            		// TODO: also need case for only idxStart (5:) or only idxEnd (:8) being specified
	                                                String[] rowRangeParts = indexString.split(":");
	                                                idxStart = Integer.valueOf(rowRangeParts[0].trim());
	                                                idxEnd   = Integer.valueOf(rowRangeParts[1].trim());
                                            	}

                                                if (idxStart >= idxEnd) {
                                                    System.out.println("error: start index >= end index in range specification for "+name+": " + indexString);
                                                } else if (
                                                		allowedArrayIndex(dim0min_val, value.length-1+dim0min_val, idxStart, name) &&
                                                		allowedArrayIndex(dim0min_val, value.length-1+dim0min_val, idxEnd,   name)) {
                                                	if (inValues.length == 1) {
                                                        // insert the same value into all possible positions
                                                    	int constantValue = Integer.valueOf(inValues[0]);
                                                        for (int idx = idxStart; idx <= idxEnd; ++idx) {
                                                            // take potentially negative minimum indices into account...
                                                            value[idx-dim0min_val] = constantValue;
                                                        }
                                                	} else {
                                                        for (int idx = idxStart; idx <= idxEnd; ++idx) {
                                                            // take potentially negative minimum indices into account...
                                                            value[idx-dim0min_val] = Integer.valueOf(inValues[idx - idxStart]);
                                                        }
                                                	}
                                                }
                                            } else {
                                                // single entry specified
                                                int idx = Integer.valueOf(indexString.trim());

                                                if (allowedArrayIndex(dim0min_val, value.length-1+dim0min_val, idx, name)) {
                                                    // take potentially negative minimum indices into account...
                                                    value[idx-dim0min_val] = Integer.valueOf(val);
                                                }
                                            }
                                        }
                                    } else {
                                        // b) whole array: myarr = 1.0, 2.0, 3.0, .... => unflatten using default (maximum) dimensions

                                        int[] flattened_value = null;
                                        if (multiIndex) {
                                            // split array values by spaces and/or commas
                                            String[] inValues = val.split("[\\s,]+");

                                            int numValues = inValues.length;
                                            List<Integer> tempValues = new LinkedList<>();

                                            for (int idx=0; idx<numValues; ++idx) {
                                                if (inValues[idx].contains("*")) {
                                                    int mi=Integer.valueOf(inValues[idx].split("\\*")[0]);
                                                    int repeatedValue = Integer.valueOf(inValues[idx].split("\\*")[1]);
                                                    for (int cnt=0; cnt<mi; ++cnt) {
                                                        tempValues.add(repeatedValue);
                                                    }
                                                } else {
                                                    tempValues.add(Integer.valueOf(inValues[idx]));
                                                }
                                            }
                                            flattened_value = tempValues.stream().mapToInt(myint->myint).toArray();
                                        } else {
                                            String[] arrayValues = val.split("[\\s,]+");
                                            flattened_value = new int[arrayValues.length];
                                            for (int idx=0; idx<arrayValues.length; ++idx) {
                                                flattened_value[idx] = Integer.valueOf(arrayValues[idx].trim());
                                            }
                                        }

                                        // unflatten until no more values leftover in Fortran order (col-first...)
                                        // assume rectangular matrix herein
                                        if (flattened_value.length <= value.length) {
                                            boolean atEnd = false;
                                            for (int idx=0; idx<value.length && !atEnd; ++idx) {
                                                if (idx<flattened_value.length) {
                                                    value[idx] = flattened_value[idx];
                                                } else {
                                                    atEnd = true;
                                                }
                                            }
                                        } else {
                                            System.out.println("ERROR: too many values specified for "+name+": " + flattened_value.length + " instead of " + value.length);
                                        }
                                    }

                                    // set field to value from namelist
                                    field.set(parseInto, value);

                                } else if (type.equals("1d:double")) {
                                    int dim0min_val = dim0min.get(name);

                                    // get field to value from namelist
                                    Field field = parseInto.getClass().getDeclaredField(fieldName);
                                    field.setAccessible(true);
                                    double[] value = (double[])field.get(parseInto);

                                    // two possibilities:
                                    if (arrayIndexSpecifierFound) {
                                        // a1) indexed value: myarr(  3) = 5.64e4  => insert at given position
                                        // a2)   range index: myarr(3:5) = 0.0     => insert into all given positions
                                    	// a3)      wildcard: myarr( : ) = 1.0, 2.0 => insert starting from beginning what is given
                                        if (indexString.contains(",")) {
                                            if (_debug) System.out.println("error: try to put value into 1d "+name+" at location indicated by more than one index: " + indexString);
                                        } else {
                                            // check for range specifiers in indices
                                            if (indexString.contains(":")) {
                                        		String[] inValues = val.split("[\\s,]+");
                                            	int idxStart;
                                            	int idxEnd;
                                            	if (indexString.length() == 1) {
                                            		// case a3) wildcard: "myarr(:) = 1.0, 2.0, 3.0"
                                            		idxStart = dim0min_val;
                                            		idxEnd = inValues.length-1 + dim0min_val;
                                            	} else {
                                            		// TODO: also need case for only idxStart (5:) or only idxEnd (:8) being specified
	                                                String[] rowRangeParts = indexString.split(":");
	                                                idxStart = Integer.valueOf(rowRangeParts[0].trim());
	                                                idxEnd   = Integer.valueOf(rowRangeParts[1].trim());
                                            	}
                                                if (idxStart >= idxEnd) {
                                                    System.out.println("error: start index >= end index in range specification for "+name+": " + indexString);
                                                } else if (
                                                		allowedArrayIndex(dim0min_val, value.length-1+dim0min_val, idxStart, name) &&
                                                		allowedArrayIndex(dim0min_val, value.length-1+dim0min_val, idxEnd,   name)) {
                                                	if (inValues.length == 1) {
                                                        // insert the same value into all possible positions
                                                    	double constantValue = Double.valueOf(inValues[0]);
                                                        for (int idx = idxStart; idx <= idxEnd; ++idx) {
                                                            // take potentially negative minimum indices into account...
                                                            value[idx-dim0min_val] = constantValue;
                                                        }
                                                	} else {
                                                        for (int idx = idxStart; idx <= idxEnd; ++idx) {
                                                            // take potentially negative minimum indices into account...
                                                            value[idx-dim0min_val] = Double.valueOf(inValues[idx - idxStart]);
                                                        }
                                                	}
                                                }
                                            } else {
                                                // single entry specified
                                                int idx = Integer.valueOf(indexString.trim());

                                                if (allowedArrayIndex(dim0min_val, value.length-1+dim0min_val, idx, name)) {
                                                    // take potentially negative minimum indices into account...
                                                    value[idx-dim0min_val] = Double.valueOf(val);
                                                }
                                            }
                                        }
                                    } else {
                                        // b) whole array: myarr = 1.0, 2.0, 3.0, .... => unflatten using default (maximum) dimensions

                                        double[] flattened_value = null;
                                        if (multiIndex) {
                                            // split array values by spaces and/or commas
                                            String[] inValues = val.split("[\\s,]+");

                                            int numValues = inValues.length;
                                            List<Double> tempValues = new LinkedList<>();

                                            for (int idx=0; idx<numValues; ++idx) {
                                                if (inValues[idx].contains("*")) {
                                                    if (_debug) System.out.println("INFO: found multi-index value: " + inValues[idx]);
                                                    int mi=Integer.valueOf(inValues[idx].split("\\*")[0]);
                                                    double repeatedValue = Double.valueOf(inValues[idx].split("\\*")[1]);
                                                    for (int cnt=0; cnt<mi; ++cnt) {
                                                        tempValues.add(repeatedValue);
                                                    }
                                                } else {
                                                    tempValues.add(Double.valueOf(inValues[idx]));
                                                }
                                            }
                                            flattened_value = tempValues.stream().mapToDouble(mydbl->mydbl).toArray();
                                        } else {
                                            String[] arrayValues = val.split("[\\s,]+");
                                            flattened_value = new double[arrayValues.length];
                                            for (int idx=0; idx<arrayValues.length; ++idx) {
                                                flattened_value[idx] = Double.valueOf(arrayValues[idx].trim());
                                            }
                                        }

                                        if (flattened_value.length <= value.length) {
                                            // unflatten until no more values leftover in Fortran order (col-first...)
                                            // assume rectangular matrix herein
                                            boolean atEnd = false;
                                            for (int idx=0; idx<value.length && !atEnd; ++idx) {
                                                if (idx<flattened_value.length) {
                                                    value[idx] = flattened_value[idx];
                                                } else {
                                                    atEnd = true;
                                                }
                                            }
                                        } else {
                                            System.out.println("ERROR: too many values specified for "+name+": " + flattened_value.length + " instead of " + value.length);
                                        }

                                    }

                                    // set field to value from namelist
                                    field.set(parseInto, value);

                                } else if (type.equals("2d:int")) {
                                    int dim0min_val = dim0min.get(name);
                                    int dim1min_val = dim1min.get(name);

                                    // get field to value from namelist
                                    Field field = parseInto.getClass().getDeclaredField(fieldName);
                                    field.setAccessible(true);
                                    int[][] value = (int[][])field.get(parseInto);

                                    // two possibilities:
                                    if (arrayIndexSpecifierFound) {
                                        // a1) indexed value: myarr(  3, 4) = 5.64e4  => insert at given position
                                        // a2)   range index: myarr(3:5, 1) = 0.0     => insert into all given positions
                                        String[] strIndices = indexString.split("[\\s]*,");
                                        if (strIndices.length != 2) {
                                            if (_debug) System.out.println("error: try to put value into 2d "+name+" at location indicated by more than two indices: " + indexString);
                                        } else {
                                            // check for range specifiers in indices
                                            if (strIndices[0].contains(":") || strIndices[1].contains(":")) {
                                                int row_start = 0, row_end = 0, col_start = 0, col_end = 0;
                                                if (strIndices[0].contains(":")) {
                                                    String[] rowRangeParts = strIndices[0].split(":");
                                                    row_start = Integer.valueOf(rowRangeParts[0].trim());
                                                    row_end   = Integer.valueOf(rowRangeParts[1].trim());
                                                } else {
                                                    row_start = Integer.valueOf(strIndices[0].trim());
                                                    row_end   = Integer.valueOf(strIndices[0].trim());
                                                }
                                                if (strIndices[1].contains(":")) {
                                                    String[] colRangeParts = strIndices[1].split(":");
                                                    col_start = Integer.valueOf(colRangeParts[0].trim());
                                                    col_end   = Integer.valueOf(colRangeParts[1].trim());
                                                } else {
                                                    col_start = Integer.valueOf(strIndices[1].trim());
                                                    col_end   = Integer.valueOf(strIndices[1].trim());
                                                }
                                                if (row_start >= row_end || col_start >= col_end) {
                                                    System.out.println("error: start index >= end index in range specification for "+name+": " + indexString);
                                                } else if (allowedArrayIndex(dim0min_val, value.length-1+dim0min_val, row_start, name)
                                                    && allowedArrayIndex(dim0min_val, value.length-1+dim0min_val, row_end, name)
                                                    && allowedArrayIndex(dim1min_val, value[0].length-1+dim1min_val, col_start, name)
                                                    && allowedArrayIndex(dim1min_val, value[0].length-1+dim1min_val, col_end, name)) {
                                                    // insert into all possible positions
                                                    for (int row = row_start; row <= row_end; ++row) {
                                                        for (int col = col_start; col <= col_end; ++col) {

                                                            // take potentially negative minimum indices into account...
                                                            value[row-dim0min_val][col-dim1min_val] = Integer.valueOf(val);
                                                        }
                                                    }
                                                }
                                            } else {
                                                // single entry specified (if single value)
                                            	// or
                                            	// starting position specified (if multiple values)

                                            	int idx_row = Integer.valueOf(strIndices[0].trim());
                                                int idx_col = Integer.valueOf(strIndices[1].trim());

                                                if (allowedArrayIndex(dim0min_val, value.length-1+dim0min_val, idx_row, name)
                                                    && allowedArrayIndex(dim1min_val, value[0].length-1+dim1min_val, idx_col, name)) {

                                                	if (val.contains(",") || val.contains(" ")) {
                                                		// multiple entries
                                                		String[] values = val.split("[\\s]*,");

                                                		for (int idxValue = 0; idxValue < values.length; ++idxValue) {

                                                			int slowOffset = idxValue / value.length;
                                                			int fastOffset = idxValue % value.length;

                                                			value[idx_row-dim0min_val + fastOffset][idx_col-dim1min_val + slowOffset] = Integer.valueOf(values[idxValue].trim());
                                                		}

                                                	} else {
	                                                    // take potentially negative minimum indices into account...
	                                                	value[idx_row-dim0min_val][idx_col-dim1min_val] = Integer.valueOf(val);
                                                	}
                                                }
                                            }
                                        }
                                    } else {
                                        // b) whole array: myarr = 1.0, 2.0, 3.0, .... => unflatten using default (maximum) dimensions

                                        int[] flattened_value = null;
                                        if (multiIndex) {
                                            // split array values by spaces and/or commas
                                            String[] inValues = val.split("[\\s,]+");

                                            int numValues = inValues.length;
                                            List<Integer> tempValues = new LinkedList<>();

                                            for (int idx=0; idx<numValues; ++idx) {
                                                if (inValues[idx].contains("*")) {
                                                    if (_debug) System.out.println("INFO: found multi-index value: " + inValues[idx]);
                                                    int mi=Integer.valueOf(inValues[idx].split("\\*")[0]);
                                                    int repeatedValue = Integer.valueOf(inValues[idx].split("\\*")[1]);
                                                    for (int cnt=0; cnt<mi; ++cnt) {
                                                        tempValues.add(repeatedValue);
                                                    }
                                                } else {
                                                    tempValues.add(Integer.valueOf(inValues[idx]));
                                                }
                                            }
                                            flattened_value = tempValues.stream().mapToInt(myint->myint).toArray();
                                        } else {
                                            String[] arrayValues = val.split("[\\s,]+");
                                            flattened_value = new int[arrayValues.length];
                                            for (int idx=0; idx<arrayValues.length; ++idx) {
                                                flattened_value[idx] = Integer.valueOf(arrayValues[idx].trim());
                                            }
                                        }

                                        if (flattened_value.length <= value.length*value[0].length) {
                                            // unflatten until no more values leftover in Fortran order (col-first...)
                                            // assume rectangular matrix herein
                                            int flatIdx = 0;
                                            boolean atEnd = false;
                                            for     (int col=0; col<value[0].length && !atEnd; ++col) {
                                                for (int row=0; row<value   .length && !atEnd; ++row) {
                                                    if (flatIdx<flattened_value.length) {
                                                        value[row][col] = flattened_value[flatIdx++];
                                                    } else {
                                                        atEnd = true;
                                                    }
                                                }
                                            }
                                        } else {
                                            System.out.println("ERROR: too many values specified for "+name+": " + flattened_value.length + " instead of " + (value.length*value[0].length));
                                        }
                                    }

                                    // set field to value from namelist
                                    field.set(parseInto, value);

                                } else if (type.equals("2d:double")) {
                                    int dim0min_val = dim0min.get(name);
                                    int dim1min_val = dim1min.get(name);

                                    // get field to value from namelist
                                    Field field = parseInto.getClass().getDeclaredField(fieldName);
                                    field.setAccessible(true);
                                    double[][] value = (double[][])field.get(parseInto);

                                    // two possibilities:
                                    if (arrayIndexSpecifierFound) {
                                        // a1) indexed value: myarr(  3, 4) = 5.64e4  => insert at given position
                                        // a2)   range index: myarr(3:5, 1) = 0.0     => insert into all given positions
                                        String[] strIndices = indexString.split("[\\s]*,");
                                        if (strIndices.length != 2) {
                                            if (_debug) System.out.println("error: try to put value into 2d "+name+" at location indicated by more than two indices: " + indexString);
                                        } else {
                                            // check for range specifiers in indices
                                            if (strIndices[0].contains(":") ||strIndices[1].contains(":")) {
                                                int row_start = 0, row_end = 0, col_start = 0, col_end = 0;
                                                if (strIndices[0].contains(":")) {
                                                    String[] rowRangeParts = strIndices[0].split(":");
                                                    row_start = Integer.valueOf(rowRangeParts[0].trim());
                                                    row_end   = Integer.valueOf(rowRangeParts[1].trim());
                                                } else {
                                                    row_start = Integer.valueOf(strIndices[0].trim());
                                                    row_end   = Integer.valueOf(strIndices[0].trim());
                                                }
                                                if (strIndices[1].contains(":")) {
                                                    String[] colRangeParts = strIndices[1].split(":");
                                                    col_start = Integer.valueOf(colRangeParts[0].trim());
                                                    col_end   = Integer.valueOf(colRangeParts[1].trim());
                                                } else {
                                                    col_start = Integer.valueOf(strIndices[1].trim());
                                                    col_end   = Integer.valueOf(strIndices[1].trim());
                                                }
                                                if (row_start >= row_end || col_start >= col_end) {
                                                    System.out.println("error: start index >= end index in range specification for "+name+": " + indexString);
                                                } else if (allowedArrayIndex(dim0min_val, value.length-1+dim0min_val, row_start, name)
                                                    && allowedArrayIndex(dim0min_val, value.length-1+dim0min_val, row_end, name)
                                                    && allowedArrayIndex(dim1min_val, value[0].length-1+dim1min_val, col_start, name)
                                                    && allowedArrayIndex(dim1min_val, value[0].length-1+dim1min_val, col_end, name)) {
                                                    // insert into all possible positions
                                                    for (int row = row_start; row <= row_end; ++row) {
                                                        for (int col = col_start; col <= col_end; ++col) {

                                                            // take potentially negative minimum indices into account...
                                                            value[row-dim0min_val][col-dim1min_val] = Double.valueOf(val);
                                                        }
                                                    }
                                                }
                                            } else {
                                                // single entry specified
                                            	// or
                                            	// starting position specified (if multiple values)

                                            	int idx_row = Integer.valueOf(strIndices[0].trim());
                                                int idx_col = Integer.valueOf(strIndices[1].trim());

                                                if (allowedArrayIndex(dim0min_val, value.length-1+dim0min_val, idx_row, name)
                                                    && allowedArrayIndex(dim1min_val, value[0].length-1+dim1min_val, idx_col, name)) {

                                                	if (val.contains(",") || val.contains(" ")) {
                                                		// multiple entries
                                                		String[] values = val.split("[\\s]*,");

                                                		for (int idxValue = 0; idxValue < values.length; ++idxValue) {

                                                			int slowOffset = idxValue / value.length;
                                                			int fastOffset = idxValue % value.length;

                                                			value[idx_row-dim0min_val + fastOffset][idx_col-dim1min_val + slowOffset] = Double.valueOf(values[idxValue].trim());
                                                		}

                                                	} else {
	                                                    // take potentially negative minimum indices into account...
	                                                	value[idx_row-dim0min_val][idx_col-dim1min_val] = Double.valueOf(val);
                                                	}
                                                }
                                            }
                                        }
                                    } else {
                                        // b) whole array: myarr = 1.0, 2.0, 3.0, .... => unflatten using default (maximum) dimensions

                                        double[] flattened_value = null;
                                        if (multiIndex) {
                                            // split array values by spaces and/or commas
                                            String[] inValues = val.split("[\\s,]+");

                                            int numValues = inValues.length;
                                            List<Double> tempValues = new LinkedList<>();

                                            for (int idx=0; idx<numValues; ++idx) {
                                                if (inValues[idx].contains("*")) {
                                                    if (_debug) System.out.println("INFO: found multi-index value: " + inValues[idx]);
                                                    int mi=Integer.valueOf(inValues[idx].split("\\*")[0]);
                                                    double repeatedValue = Double.valueOf(inValues[idx].split("\\*")[1]);
                                                    for (int cnt=0; cnt<mi; ++cnt) {
                                                        tempValues.add(repeatedValue);
                                                    }
                                                } else {
                                                    tempValues.add(Double.valueOf(inValues[idx]));
                                                }
                                            }
                                            flattened_value = tempValues.stream().mapToDouble(mydbl->mydbl).toArray();
                                        } else {
                                            String[] arrayValues = val.split("[\\s,]+");
                                            flattened_value = new double[arrayValues.length];
                                            for (int idx=0; idx<arrayValues.length; ++idx) {
                                                flattened_value[idx] = Double.valueOf(arrayValues[idx].trim());
                                            }
                                        }

                                        if (flattened_value.length <= value.length*value[0].length) {
                                            // unflatten until no more values leftover in Fortran order (col-first...)
                                            // assume rectangular matrix herein
                                            int flatIdx = 0;
                                            boolean atEnd = false;
                                            for     (int col=0; col<value[0].length && !atEnd; ++col) {
                                                for (int row=0; row<value   .length && !atEnd; ++row) {
                                                    if (flatIdx<flattened_value.length) {
                                                        value[row][col] = flattened_value[flatIdx++];
                                                    } else {
                                                        atEnd = true;
                                                    }
                                                }
                                            }
                                        } else {
                                            System.out.println("ERROR: too many values specified for "+name+": " + flattened_value.length + " instead of " + (value.length*value[0].length));
                                        }
                                    }

                                    // set field to value from namelist
                                    field.set(parseInto, value);

                                } else {
                                    throw new RuntimeException("Type " + type + " is not implemented yet in FortranNamelist parser!");
                                }
                            } catch (Exception e) {
                                throw new RuntimeException(e);
                            }
                        }

                        // next variable name
                        if (end<part.length()) {
                            name = part.substring(end+1).trim().toLowerCase();
                            // some namelists have their entries separated by ",", so remove the commas in this step
                            if (name.startsWith(",")) {
                                name = name.substring(1).trim();
                            }
                            if (_debug) System.out.println("next name: " + name);
                        }
                    }

                    // look for forgotten variable values: leave at default values
                    else if (names.containsKey(part.toLowerCase())) {
                        System.out.println("INFO: probably forgot value for \"" + name + "\"? [only found next variable name \""+part.toLowerCase()+"\"] ==> will use default values");
                    }

                    // You should not see this, since it indicates that a certain part could not be interpreted.
                    else {
                        System.out.println("last name: " + name.toLowerCase());
                        System.out.println("part from "+last_i+ " to " + i + " : \"" + part+"\"");
                    }

                }
                i++;

                var_counter++;
            }
        } else {
        	throw new RuntimeException("namelist '"+groupName+"' not found in the given input");
        }
    }

    /**
     * Check index boundaries for inserting values into array; named version.
     *
     * @param minIdx     minimum allowable index, i.e. 0 for int[] a = int[10];
     * @param maxIdx     maximum allowable index, i.e. 9 for int[] a = int[10];
     * @param idxToCheck index to check, i.e. 5 or 11
     * @param varname    name of variable to be included in error message for
     *                   debugging
     * @return true: idxToCheck within bounds (i.e. for idxToCheck=5); false if not
     *         (i.e. for idxToCheck=11)
     */
    public static boolean allowedArrayIndex(int minIdx, int maxIdx, int idxToCheck, String varname) {
        if (idxToCheck < minIdx || idxToCheck > maxIdx) {
            System.out.println("ERROR: try to insert value into " + varname + " at " + idxToCheck + ", allowed: ["
                + minIdx + ":" + maxIdx + "]");
            return false;
        }
        return true;
    }

    // ---------------------

    public static int[] stripTrailingZeros(int[] input) {
		if (input == null) {
			return null;
		}

		int nonZeroLength = input.length;
		if (nonZeroLength == 0) {
			return new int[0];
		}

		while (nonZeroLength > 0 && input[nonZeroLength - 1] == 0) {
			nonZeroLength--;
		}
		return Arrays.copyOf(input, nonZeroLength);
	}

	public static double[] stripTrailingZeros(double[] input) {
		final int minLength = 0;
		return stripTrailingZeros(input, minLength);
	}

	public static double[] stripTrailingZeros(double[] input, int minLength) {
		if (input == null) {
			return null;
		}

		int nonZeroLength = input.length;
		if (nonZeroLength == 0) {
			return new double[0];
		}

		while (nonZeroLength > 0 && input[nonZeroLength - 1] == 0) {
			nonZeroLength--;
		}

		if (nonZeroLength > 0 || minLength == 0) {
			return Arrays.copyOf(input, nonZeroLength);
		} else {
			return new double[minLength];
		}
	}

	// ---------------------

	public static final String toFortranAsString(final boolean b) {
		return b ? ".true." : ".false.";
	}

	public static final String toFortranAsString(final double[] x) {
		final boolean omitZerosAtEnd = true;
		return toFortranAsString(x, omitZerosAtEnd);
	}

	public static final String toFortranAsString(final double[] x, final boolean omitZerosAtEnd) {
		Objects.requireNonNull(x);
		String s = "";
		for (int i = 0; i < x.length; ++i) {
			s += String.format(Locale.ENGLISH, "% .20e", x[i]);

			boolean onlyZerosFollow = true;
			if (omitZerosAtEnd) {
				for (int j = i + 1; j < x.length; ++j) {
					if (x[j] != 0.0) {
						onlyZerosFollow = false;
						break;
					}
				}
			}

			if (i < x.length - 1) {
				if (onlyZerosFollow) {
					if (omitZerosAtEnd) {
						break;
					} else {
						s += ", ";
					}
				} else {
					s += ", ";
				}
			} // no matter what, the last number is not followed by a comma
		}

		return s;
	}

	public static final String toFortranAsString(final int[] x) {
		final boolean omitZerosAtEnd = true;
		return toFortranAsString(x, omitZerosAtEnd);
	}

	public static final String toFortranAsString(final int[] x, final boolean omitZerosAtEnd) {
		Objects.requireNonNull(x);
		String s = "";
		for (int i = 0; i < x.length; ++i) {
			s += String.format(Locale.ENGLISH, "%d", x[i]);

			boolean onlyZerosFollow = true;
			if (omitZerosAtEnd) {
				for (int j = i + 1; j < x.length; ++j) {
					if (x[j] != 0) {
						onlyZerosFollow = false;
						break;
					}
				}
			}

			if (i < x.length - 1) {
				if (onlyZerosFollow) {
					if (omitZerosAtEnd) {
						break;
					} else {
						s += ", ";
					}
				} else {
					s += ", ";
				}
			} // no matter what, the last number is not followed by a comma
		}

		return s;
	}
}