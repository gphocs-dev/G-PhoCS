/*
 * This class' methods are used for validation purposes
 */
package CFG.ControlFileGenerator;

/**
 *
 * @author Tal
 */
public class Validate {
    
    // Validates that a given input string is a positive integer
    public static boolean validateInt(String input) {
        if("".equals(input))
            return true;
        int len = input.length();
        for (int i = 0; i < len; i++) {
            if(!Character.isDigit(input.charAt(i))) {
                return false;
            }
        }
        return true;
    }

    // Validates that a given input string is a positive numerical number
    public static boolean validateDouble(String input) {
        if("".equals(input))
            return true;
        if(".".equals(input))
            return false;
        int len = input.length();
        int dotCounter = 0;
        for (int i = 0; i < len; i++) {
            if(!Character.isDigit(input.charAt(i))) {
                if (input.charAt(i) != '.') {
                    return false;
                } else {
                    if (dotCounter == 0) {
                        dotCounter++;
                    }
                    else {
                        return false;
                    }
                }
            }
        }
        return true;
    }
    
    // Helper function to check invalid input for the tree.
    public static boolean validateTreeInput(String treeInput) {
        int len = treeInput.length();
        if (len == 0) {
            return false;
        }
        int counterOpen = 0;
        int counterClose = 0;
        int counterComma = 0;
        char curSign;
        char lastSign = treeInput.charAt(len - 1);
        if (lastSign == '(' || lastSign == ')' || lastSign == ',') {
            // Checks if the last sign is not a braket or a comma
            return true;
        }
        for (int i = 0; i < len; i++) {
            curSign = treeInput.charAt(i);
            if (curSign == '(') {
                counterOpen++;
            } else if (curSign == ')') {
                counterClose++;
            } else if (curSign == ',') {
                counterComma++;
            }

            // Checks if there are more closing brackets before open ones
            if (counterClose > counterOpen) {
                return true;
            }
        }
        // Checks that the number of open brackets and closing brackets is equal
        if (counterOpen != counterClose) {
            return true;
            // Checks that the number of commas is equal to the number of brackets
        } else if (counterOpen != counterComma) {
            return true;
        }

        return false;
    }
    
    public static boolean validateNewickTree(String treeInput, int inputLen) {
        int commaCounter = 0;
        NewickTree nwtValidate = new NewickTree(treeInput);
        
        // Checks there is no whitespace in one of the names
        for (int i = 0; i < nwtValidate.childArrayLen; i++) {
            if (checkStringForWhiteSpace(nwtValidate.childArray[i].data))
                return true;
        }
        
        for (int i = 0; i < nwtValidate.parentArrayLen; i++) {
            if (checkStringForWhiteSpace(nwtValidate.parentArray[i].data))
                return true;
        }

        
        // Checks for duplicates
        for (int i = 0; i < nwtValidate.childArrayLen; i++) {
            for (int j = 0; j < nwtValidate.childArrayLen; j++) {
                if (nwtValidate.childArray[i].data.compareTo(nwtValidate.childArray[j].data) == 0 && i != j)
                    return true;
            }
            for (int j = 0; j < nwtValidate.parentArrayLen; j++) {
                if (nwtValidate.childArray[i].data.compareTo(nwtValidate.parentArray[j].data) == 0)
                    return true;
            }
        }
        
        for (int i = 0; i < nwtValidate.parentArrayLen; i++) {
            for (int j = 0; j < nwtValidate.childArrayLen; j++) {
                if (nwtValidate.parentArray[i].data.compareTo(nwtValidate.childArray[j].data) == 0)
                    return true;
            }
            for (int j = 0; j < nwtValidate.parentArrayLen; j++) {
                if (nwtValidate.parentArray[i].data.compareTo(nwtValidate.parentArray[j].data) == 0 && i != j)
                    return true;
            }
        }
        
        for (int i = 0; i < treeInput.length(); i ++) {
            char c = treeInput.charAt(i);
            if(c == ',')
                commaCounter++;
        }
        
        // Checking if the number of commas fits the given input
        if (commaCounter * 2 + 1 != nwtValidate.childArrayLen + nwtValidate.parentArrayLen)
            return true;
        
        return false;
    }
    
    private static boolean checkStringForWhiteSpace(String s) {
        for(int i = 0; i < s.length(); i++){
            if(Character.isWhitespace(s.charAt(i))){
                return true;
            }
        }
        return false;
    }
}
