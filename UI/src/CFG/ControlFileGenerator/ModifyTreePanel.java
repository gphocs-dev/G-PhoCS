/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package CFG.ControlFileGenerator;

import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.SpringLayout;

/**
 *
 * @author Tal
 */
public class ModifyTreePanel {
    
    public static void modifyCurrentLabel(BSTNode[] nodeArray, String[] inputArray) {
        
        int numOfFields = nodeArray.length;
        JTextField[] textFields = new JTextField[numOfFields];
        String[] labels = new String[numOfFields];
        for (int i = 0; i < numOfFields; i++) {
            textFields[i] = new JTextField(5);
            textFields[i].setText(inputArray[i]);
            labels[i] = "Insert the sample of Current POP " + nodeArray[i].data + ":";
        }
        
        JPanel inputPanel = new JPanel(new SpringLayout());
        for (int i = 0; i < numOfFields; i++) {
            JLabel label = new JLabel(labels[i], JLabel.TRAILING);
            inputPanel.add(label);
            label.setLabelFor(textFields[i]);
            inputPanel.add(textFields[i]);
        }
        
        SpringUtilities.makeCompactGrid(inputPanel,
                                numOfFields, 2, //rows, cols
                                6, 6,        //initX, initY
                                6, 6);       //xPad, yPad
        int result = JOptionPane.showConfirmDialog(null, inputPanel,"Insert the Samples of the Current-POP", JOptionPane.OK_CANCEL_OPTION);
        if (result == JOptionPane.OK_OPTION) {
            for (int i = 0; i < textFields.length; i++) {
                inputArray[i] = textFields[i].getText();
            }
        }
        
    }
    
        public static void modifyAncestralTauInitial(BSTNode[] nodeArray, String[] inputArray) {
        
        int numOfFields = nodeArray.length;
        JTextField[] textFields = new JTextField[numOfFields];
        String[] labels = new String[numOfFields];
        String[] childrenLabels = new String[numOfFields];
        for (int i = 0; i < numOfFields; i++) {
            textFields[i] = new JTextField(5);
            textFields[i].setText(inputArray[i]);
            labels[i] = "Insert the Tau-Initial of Ancestral-POP " + nodeArray[i].data + ":";
//            childrenLabels[i] = "Daughters of Ancestral-POP: " + nodeArray[i].left.data + " , " + nodeArray[i].right.data;
        }
        
        JPanel inputPanel = new JPanel(new SpringLayout());
        for (int i = 0; i < numOfFields; i++) {
            JLabel label = new JLabel(labels[i], JLabel.TRAILING);
            inputPanel.add(label);
            label.setLabelFor(textFields[i]);
            inputPanel.add(textFields[i]);
        }
        
        SpringUtilities.makeCompactGrid(inputPanel,
                                numOfFields, 2, //rows, cols
                                6, 6,        //initX, initY
                                6, 6);       //xPad, yPad
        int result = JOptionPane.showConfirmDialog(null, inputPanel,"Insert the Tau-Initial of the Ancestral-POPs", JOptionPane.OK_CANCEL_OPTION);
        if (result == JOptionPane.OK_OPTION) {
            for (int i = 0; i < textFields.length; i++) {
                inputArray[i] = textFields[i].getText();
            }
        }
        
    }
    
}
