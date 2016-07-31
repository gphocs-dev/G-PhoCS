/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package CFG.ControlFileGenerator;

import java.awt.Color;
import javax.swing.ButtonGroup;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTextField;
import javax.swing.SpringLayout;

/**
 *
 * @author Tal
 */
public class InputPanel {

    public static String[] POP_START = {"POP-START", "", "", "", ""};
    public static String[] POP_END = {"POP-END", "", "", "", ""};
    public static String NAME = "name";
    public static String SAMPLE = "samples";
    public static String CHILDREN = "children";
    public static String TAU_INITIAL = "tau-initial";
    public static String HAPLOID = "haploid";
    public static String DIPLOID = "diploid";

    public static boolean insertCurrentLabel(BSTNode[] nodeArray, String[] inputArray, String[] hd) {

        boolean operationSuccess = false;
        int numOfFields = nodeArray.length;
        JTextField[] nameTextFields = new JTextField[numOfFields];
        JTextField[] labelTextFields = new JTextField[numOfFields];
        ButtonGroup[] haploidDiploidGroup = new ButtonGroup[numOfFields];
        JRadioButton[] haploid = new JRadioButton[numOfFields];
        JRadioButton[] diploid = new JRadioButton[numOfFields];

        for (int i = 0; i < numOfFields; i++) {
            nameTextFields[i] = new JTextField(5);
            nameTextFields[i].setEnabled(false);
            nameTextFields[i].setText(nodeArray[i].data);
            nameTextFields[i].setDisabledTextColor(Color.black);
            labelTextFields[i] = new JTextField(5);
            labelTextFields[i].setText(inputArray[i]);
            haploidDiploidGroup[i] = new ButtonGroup();
            haploid[i] = new JRadioButton();
            diploid[i] = new JRadioButton();
            haploidDiploidGroup[i].add(haploid[i]);
            haploidDiploidGroup[i].add(diploid[i]);
            haploid[i].setText(HAPLOID);
            diploid[i].setText(DIPLOID);
        }

        JPanel inputPanel = new JPanel(new SpringLayout());
        for (int i = 0; i < numOfFields; i++) {
            for (int j = 0; j < 5; j++) {
                JLabel label = new JLabel(POP_START[j], JLabel.TRAILING);
                inputPanel.add(label);
            }

            JLabel label = new JLabel("", JLabel.TRAILING);
            inputPanel.add(label);
            label = new JLabel(NAME, JLabel.TRAILING);
            inputPanel.add(label);
            label.setLabelFor(nameTextFields[i]);
            inputPanel.add(nameTextFields[i]);
            label = new JLabel("", JLabel.TRAILING);
            inputPanel.add(label);
            label = new JLabel("", JLabel.TRAILING);
            inputPanel.add(label);

            label = new JLabel("", JLabel.TRAILING);
            inputPanel.add(label);
            label = new JLabel(SAMPLE, JLabel.TRAILING);
            inputPanel.add(label);
            label.setLabelFor(labelTextFields[i]);
            inputPanel.add(labelTextFields[i]);
            label.setLabelFor(haploid[i]);
            inputPanel.add(haploid[i]);
            label.setLabelFor(diploid[i]);
            inputPanel.add(diploid[i]);
            if ("h".equals(hd[i])) {
                haploid[i].setSelected(true);
            } else {
                diploid[i].setSelected(true);
            }

            for (int j = 0; j < 5; j++) {
                label = new JLabel(POP_END[j], JLabel.TRAILING);
                inputPanel.add(label);
            }
        }

        SpringUtilities.makeCompactGrid(inputPanel,
                numOfFields * 4, 5, //rows, cols
                6, 6, //initX, initY
                6, 6);       //xPad, yPad
        int result = JOptionPane.showConfirmDialog(null, inputPanel, "Insert the Samples of the Current-POP", JOptionPane.OK_CANCEL_OPTION);

        boolean noEmptySamples = true;
        for (int i = 0; i < labelTextFields.length; i++) {
            if("".equals(labelTextFields[i].getText()))
            {
                noEmptySamples = false;
                break;
            }
        }
        
        boolean noIdenticalSamples = true;
        for (int i = 0; i < labelTextFields.length; i++) {
            for (int j = 0; j < labelTextFields.length; j++) {
                if (!("".equals(labelTextFields[i].getText()))) {
                    if (i != j && labelTextFields[i].getText().equals(labelTextFields[j].getText())) {
                        noIdenticalSamples = false;
                        break;
                    }
                }
            }
        }
        if (result == JOptionPane.OK_OPTION) {
            ControlFileGeneratorUI.samplesAll.clear();
            if (!noIdenticalSamples) {
                JOptionPane.showMessageDialog(null, "Can't have two different pops with identical samples", "message", JOptionPane.ERROR_MESSAGE);
                operationSuccess = insertCurrentLabelAfterError(nodeArray, inputArray, hd, nameTextFields, labelTextFields, haploidDiploidGroup, haploid, diploid);
            } else if (!noEmptySamples) {
                JOptionPane.showMessageDialog(null, "Can't have empty samples", "message", JOptionPane.ERROR_MESSAGE);
                operationSuccess = insertCurrentLabelAfterError(nodeArray, inputArray, hd, nameTextFields, labelTextFields, haploidDiploidGroup, haploid, diploid);
            } else {
                for (int i = 0; i < labelTextFields.length; i++) {
                    inputArray[i] = labelTextFields[i].getText();
                    if (diploid[i].isSelected()) {
                        hd[i] = "d";
                    } else {
                        hd[i] = "h";
                    }
                }
                
                if(CheckAllSamples(inputArray))
                    operationSuccess = true;
                else {
                    JOptionPane.showMessageDialog(null, "Can't have two different pops with identical samples", "message", JOptionPane.ERROR_MESSAGE);
                    operationSuccess = insertCurrentLabelAfterError(nodeArray, inputArray, hd, nameTextFields, labelTextFields, haploidDiploidGroup, haploid, diploid);
                }
                    
            }
        }
        return operationSuccess;
    }
    
    
    
    public static boolean insertCurrentLabelAfterError(BSTNode[] nodeArray, String[] inputArray, String[] hd, JTextField[] nameTextFields, JTextField[] labelTextFields,
                                                            ButtonGroup[] haploidDiploidGroup, JRadioButton[] haploid, JRadioButton[] diploid) {

        boolean operationSuccess = false;
        int numOfFields = nodeArray.length;

        JPanel inputPanel = new JPanel(new SpringLayout());
        for (int i = 0; i < numOfFields; i++) {
            for (int j = 0; j < 5; j++) {
                JLabel label = new JLabel(POP_START[j], JLabel.TRAILING);
                inputPanel.add(label);
            }

            JLabel label = new JLabel("", JLabel.TRAILING);
            inputPanel.add(label);
            label = new JLabel(NAME, JLabel.TRAILING);
            inputPanel.add(label);
            label.setLabelFor(nameTextFields[i]);
            inputPanel.add(nameTextFields[i]);
            label = new JLabel("", JLabel.TRAILING);
            inputPanel.add(label);
            label = new JLabel("", JLabel.TRAILING);
            inputPanel.add(label);

            label = new JLabel("", JLabel.TRAILING);
            inputPanel.add(label);
            label = new JLabel(SAMPLE, JLabel.TRAILING);
            inputPanel.add(label);
            label.setLabelFor(labelTextFields[i]);
            inputPanel.add(labelTextFields[i]);
            label.setLabelFor(haploid[i]);
            inputPanel.add(haploid[i]);
            label.setLabelFor(diploid[i]);
            inputPanel.add(diploid[i]);
            if ("h".equals(hd[i])) {
                haploid[i].setSelected(true);
            } else {
                diploid[i].setSelected(true);
            }

            for (int j = 0; j < 5; j++) {
                label = new JLabel(POP_END[j], JLabel.TRAILING);
                inputPanel.add(label);
            }
        }

        SpringUtilities.makeCompactGrid(inputPanel,
                numOfFields * 4, 5, //rows, cols
                6, 6, //initX, initY
                6, 6);       //xPad, yPad
        int result = JOptionPane.showConfirmDialog(null, inputPanel, "Insert the Samples of the Current-POP", JOptionPane.OK_CANCEL_OPTION);

        boolean noEmptySamples = true;
        for (int i = 0; i < labelTextFields.length; i++) {
            if("".equals(labelTextFields[i].getText()))
            {
                noEmptySamples = false;
                break;
            }
        }
        
        boolean noIdenticalSamples = true;
        for (int i = 0; i < labelTextFields.length; i++) {
            for (int j = 0; j < labelTextFields.length; j++) {
                if (!("".equals(labelTextFields[i].getText()))) {
                    if (i != j && labelTextFields[i].getText().equals(labelTextFields[j].getText())) {
                        noIdenticalSamples = false;
                        break;
                    }
                }
            }
        }
        if (result == JOptionPane.OK_OPTION) {
            ControlFileGeneratorUI.samplesAll.clear();
            if (!noIdenticalSamples) {
                JOptionPane.showMessageDialog(null, "Can't have two different pops with identical samples", "message", JOptionPane.ERROR_MESSAGE);
                operationSuccess = insertCurrentLabelAfterError(nodeArray, inputArray, hd, nameTextFields, labelTextFields, haploidDiploidGroup, haploid, diploid);
            } else if (!noEmptySamples) {
                JOptionPane.showMessageDialog(null, "Can't have empty samples", "message", JOptionPane.ERROR_MESSAGE);
                operationSuccess = insertCurrentLabelAfterError(nodeArray, inputArray, hd, nameTextFields, labelTextFields, haploidDiploidGroup, haploid, diploid);
            } else {
                for (int i = 0; i < labelTextFields.length; i++) {
                    inputArray[i] = labelTextFields[i].getText();
                    if (diploid[i].isSelected()) {
                        hd[i] = "d";
                    } else {
                        hd[i] = "h";
                    }
                }
                if(CheckAllSamples(inputArray))
                    operationSuccess = true;
                else {
                    JOptionPane.showMessageDialog(null, "Can't have two different pops with identical samples", "message", JOptionPane.ERROR_MESSAGE);
                    operationSuccess = insertCurrentLabelAfterError(nodeArray, inputArray, hd, nameTextFields, labelTextFields, haploidDiploidGroup, haploid, diploid);
                }
            }
        }
        return operationSuccess;
    }
        
    public static boolean CheckAllSamples(String[] inputArray)
    {
        for(int i = 0; i < inputArray.length; i++)
        {
            String[] splited = inputArray[i].split("\\s+");
            for(int j = 0; j < splited.length; j++)
            {
                ControlFileGeneratorUI.samplesAll.add(splited[j]);
            }
        }
        
   //     return true;
        
        boolean noDuplicate = true;
        String curSample;
        String[] allSamples = ControlFileGeneratorUI.samplesAll.toArray(new String[ControlFileGeneratorUI.samplesAll.size()]);
        for(int i = 0; i < allSamples.length; i++)
        {
            curSample = allSamples[i];
            for(int j = 0; j < allSamples.length; j++)
            {
                if(j != i && curSample.equals(allSamples[j]))
                    noDuplicate = false;
            }
        }
        
        return noDuplicate;
    }
        
        

    public static boolean insertAncestralTauInitial(BSTNode[] nodeArray, String[] inputArray) {

        boolean operationSuccess = false;
        int numOfFields = nodeArray.length;
        JTextField[] nameTextFields = new JTextField[numOfFields];
        JTextField[] leftChildTextFields = new JTextField[numOfFields];
        JTextField[] rightChildTextFields = new JTextField[numOfFields];
        JTextField[] tauTextFields = new JTextField[numOfFields];
        for (int i = 0; i < numOfFields; i++) {
            nameTextFields[i] = new JTextField(5);
            nameTextFields[i].setEnabled(false);
            nameTextFields[i].setText(nodeArray[i].data);
            nameTextFields[i].setDisabledTextColor(Color.black);
            leftChildTextFields[i] = new JTextField(5);
            leftChildTextFields[i].setEnabled(false);
            leftChildTextFields[i].setText(nodeArray[i].left.data);
            leftChildTextFields[i].setDisabledTextColor(Color.black);
            rightChildTextFields[i] = new JTextField(5);
            rightChildTextFields[i].setEnabled(false);
            rightChildTextFields[i].setText(nodeArray[i].right.data);
            rightChildTextFields[i].setDisabledTextColor(Color.black);
            tauTextFields[i] = new JTextField(5);
            tauTextFields[i].setText(inputArray[i]);
        }

        JPanel inputPanel = new JPanel(new SpringLayout());
        for (int i = 0; i < numOfFields; i++) {
            for (int j = 0; j < 4; j++) {
                JLabel label = new JLabel(POP_START[j], JLabel.TRAILING);
                inputPanel.add(label);
            }

            JLabel label = new JLabel("", JLabel.TRAILING);
            inputPanel.add(label);
            label = new JLabel(NAME, JLabel.TRAILING);
            inputPanel.add(label);
            label.setLabelFor(nameTextFields[i]);
            inputPanel.add(nameTextFields[i]);
            label = new JLabel("", JLabel.TRAILING);
            inputPanel.add(label);

            label = new JLabel("", JLabel.TRAILING);
            inputPanel.add(label);
            label = new JLabel(CHILDREN, JLabel.TRAILING);
            inputPanel.add(label);
            label.setLabelFor(leftChildTextFields[i]);
            inputPanel.add(leftChildTextFields[i]);
            label.setLabelFor(rightChildTextFields[i]);
            inputPanel.add(rightChildTextFields[i]);

            label = new JLabel("", JLabel.TRAILING);
            inputPanel.add(label);
            label = new JLabel(TAU_INITIAL, JLabel.TRAILING);
            inputPanel.add(label);
            label.setLabelFor(tauTextFields[i]);
            inputPanel.add(tauTextFields[i]);
            label = new JLabel("", JLabel.TRAILING);
            inputPanel.add(label);

            for (int j = 0; j < 4; j++) {
                label = new JLabel(POP_END[j], JLabel.TRAILING);
                inputPanel.add(label);
            }
        }

        SpringUtilities.makeCompactGrid(inputPanel,
                numOfFields * 5, 4, //rows, cols
                6, 6, //initX, initY
                6, 6);       //xPad, yPad

        int result = JOptionPane.showConfirmDialog(null, inputPanel, "Insert the Tau-Initial of the Ancestral-POPs", JOptionPane.OK_CANCEL_OPTION);
        if (result == JOptionPane.OK_OPTION) {
            boolean allTauCanBeParsed = true;
            for (int i = 0; i < tauTextFields.length; i++) {
                if (!Validate.validateDouble(tauTextFields[i].getText())) {
                    allTauCanBeParsed = false;
                    break;
                }
            }
            if (!allTauCanBeParsed) {
                JOptionPane.showMessageDialog(null, "all tau-initial must be a number!", "message", JOptionPane.ERROR_MESSAGE);
                operationSuccess = insertAncestralTauInitialAfterError(nodeArray, inputArray, nameTextFields, leftChildTextFields, rightChildTextFields, tauTextFields);
            } else {
                boolean tauIsLegal = true;
                for (int i = 0; i < tauTextFields.length; i++) {
                    if (!checkIfTauSmallerThanParent(nodeArray[i], nodeArray, tauTextFields[i], tauTextFields)) {
                        tauIsLegal = false;
                        break;
                    }
                }
                if (tauIsLegal) {
                    for (int i = 0; i < tauTextFields.length; i++) {
                        inputArray[i] = tauTextFields[i].getText();
                    }
                    operationSuccess = true;
                } else {
                    JOptionPane.showMessageDialog(null, "The tau-initial of a daughter population can't be larger than its parent!", "message", JOptionPane.ERROR_MESSAGE);
                    operationSuccess = insertAncestralTauInitialAfterError(nodeArray, inputArray, nameTextFields, leftChildTextFields, rightChildTextFields, tauTextFields);
                }
            }
        }
        return operationSuccess;
    }
    
    
    
    
    public static boolean insertAncestralTauInitialAfterError(BSTNode[] nodeArray, String[] inputArray, JTextField[] nameTextFields, JTextField[] leftChildTextFields,
                                                              JTextField[] rightChildTextFields, JTextField[] tauTextFields) {

        boolean operationSuccess = false;
        int numOfFields = nodeArray.length;

        JPanel inputPanel = new JPanel(new SpringLayout());
        for (int i = 0; i < numOfFields; i++) {
            for (int j = 0; j < 4; j++) {
                JLabel label = new JLabel(POP_START[j], JLabel.TRAILING);
                inputPanel.add(label);
            }

            JLabel label = new JLabel("", JLabel.TRAILING);
            inputPanel.add(label);
            label = new JLabel(NAME, JLabel.TRAILING);
            inputPanel.add(label);
            label.setLabelFor(nameTextFields[i]);
            inputPanel.add(nameTextFields[i]);
            label = new JLabel("", JLabel.TRAILING);
            inputPanel.add(label);

            label = new JLabel("", JLabel.TRAILING);
            inputPanel.add(label);
            label = new JLabel(CHILDREN, JLabel.TRAILING);
            inputPanel.add(label);
            label.setLabelFor(leftChildTextFields[i]);
            inputPanel.add(leftChildTextFields[i]);
            label.setLabelFor(rightChildTextFields[i]);
            inputPanel.add(rightChildTextFields[i]);

            label = new JLabel("", JLabel.TRAILING);
            inputPanel.add(label);
            label = new JLabel(TAU_INITIAL, JLabel.TRAILING);
            inputPanel.add(label);
            label.setLabelFor(tauTextFields[i]);
            inputPanel.add(tauTextFields[i]);
            label = new JLabel("", JLabel.TRAILING);
            inputPanel.add(label);

            for (int j = 0; j < 4; j++) {
                label = new JLabel(POP_END[j], JLabel.TRAILING);
                inputPanel.add(label);
            }
        }

        SpringUtilities.makeCompactGrid(inputPanel,
                numOfFields * 5, 4, //rows, cols
                6, 6, //initX, initY
                6, 6);       //xPad, yPad

        int result = JOptionPane.showConfirmDialog(null, inputPanel, "Insert the Tau-Initial of the Ancestral-POPs", JOptionPane.OK_CANCEL_OPTION);
        if (result == JOptionPane.OK_OPTION) {
            boolean allTauCanBeParsed = true;
            for (int i = 0; i < tauTextFields.length; i++) {
                if (!Validate.validateDouble(tauTextFields[i].getText())) {
                    allTauCanBeParsed = false;
                    break;
                }
            }
            if (!allTauCanBeParsed) {
                JOptionPane.showMessageDialog(null, "all tau-initial must be a number!", "message", JOptionPane.ERROR_MESSAGE);
                operationSuccess = insertAncestralTauInitialAfterError(nodeArray, inputArray, nameTextFields, leftChildTextFields, rightChildTextFields, tauTextFields);
            } else {
                boolean tauIsLegal = true;
                for (int i = 0; i < tauTextFields.length; i++) {
                    if (!checkIfTauSmallerThanParent(nodeArray[i], nodeArray, tauTextFields[i], tauTextFields)) {
                        tauIsLegal = false;
                        break;
                    }
                }
                if (tauIsLegal) {
                    for (int i = 0; i < tauTextFields.length; i++) {
                        inputArray[i] = tauTextFields[i].getText();
                    }
                    operationSuccess = true;
                } else {
                    JOptionPane.showMessageDialog(null, "The tau-initial of a daughter population can't be larger than its parent!", "message", JOptionPane.ERROR_MESSAGE);
                    operationSuccess = insertAncestralTauInitialAfterError(nodeArray, inputArray, nameTextFields, leftChildTextFields, rightChildTextFields, tauTextFields);
                }
            }
        }
        return operationSuccess;
    }

    private static boolean checkIfTauSmallerThanParent(BSTNode node, BSTNode[] nodeArray, JTextField currentTauField, JTextField[] tau) {
        if (node.parent != null && !("".equals(currentTauField.getText()))) {
            for (int i = 0; i < nodeArray.length; i++) {
                if (node.parent.equals(nodeArray[i])) {
                    if (!("".equals(tau[i].getText()))) {
                        double currentTau = Double.parseDouble(currentTauField.getText());
                        double parentTau = Double.parseDouble(tau[i].getText());
                        if (currentTau >= parentTau) {
                            return false;
                        }
                    }
                }
            }
        }
        return true;
    }

}
