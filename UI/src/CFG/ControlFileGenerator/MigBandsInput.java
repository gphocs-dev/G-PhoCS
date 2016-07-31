/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package CFG.ControlFileGenerator;

import java.awt.Color;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.SpringLayout;

/**
 *
 * @author Tal
 */
public class MigBandsInput {
    
    final static String SOURCE = "source";
    final static String TARGET = "target";
    final static int NUM_OF_FIELDS = 2;

    
    public static void AddNewMigBand(BSTNode[] currentPopArray, BSTNode[] ancestralPopArray, String[] migBandNewInput, String sourceMigBand, String targetMigBand) {
        
            if (sourceMigBand.equals(targetMigBand)) {
                JOptionPane.showMessageDialog(null, "source and target can't be the same population", "message", JOptionPane.ERROR_MESSAGE);
            } else if (sourceTargetPairingExists(sourceMigBand, targetMigBand)){
                JOptionPane.showMessageDialog(null, "A migration band from " + sourceMigBand + " to " + targetMigBand + " already exists" , "message", JOptionPane.ERROR_MESSAGE);
            } else if (!popIsDaughterOfAncestral(sourceMigBand, targetMigBand, currentPopArray, ancestralPopArray)) {
                JOptionPane.showMessageDialog(null, "the chosen populations can't be a daughter and an ancestral", "message", JOptionPane.ERROR_MESSAGE); 
            } else {    
                migBandNewInput[0] = sourceMigBand;
                migBandNewInput[1] = targetMigBand;
            }
        }

    
//    public static int AddNewMigBand(BSTNode[] currentPopArray, BSTNode[] ancestralPopArray, String[] migBandNewInput, int onComplete) {
//        
//        JLabel label;
//        JTextField[] textFields = new JTextField[NUM_OF_FIELDS];
//        String[] migBandInput = {"", ""};
//
//        for (int i = 0; i < NUM_OF_FIELDS; i++) {
//            textFields[i] = new JTextField(5);
//        }
//
//        JPanel addMigBandPanel = new JPanel(new SpringLayout());
//        
//        
//        label = new JLabel(SOURCE, JLabel.TRAILING);
//        addMigBandPanel.add(label);
//        label.setLabelFor(textFields[0]);
//        addMigBandPanel.add(textFields[0]);
//        label = new JLabel(TARGET, JLabel.TRAILING);
//        addMigBandPanel.add(label);
//        label.setLabelFor(textFields[1]);
//        addMigBandPanel.add(textFields[1]);             
//        
//
//        SpringUtilities.makeCompactGrid(addMigBandPanel,
//                
//                
//                1, NUM_OF_FIELDS * 2, //rows, cols
//                6, 6, //initX, initY
//                6, 6);       //xPad, yPad
//
//        int result = JOptionPane.showConfirmDialog(null, addMigBandPanel, "Insert the source, target and mig-rate-print of the Migration-Band", JOptionPane.OK_CANCEL_OPTION);
//
//        if (result == JOptionPane.OK_OPTION) {
//            if ("".equals(textFields[0].getText()) || "".equals(textFields[1].getText())) {
//                JOptionPane.showMessageDialog(null, "One or more of the fields is empty", "message", JOptionPane.ERROR_MESSAGE);
//                onComplete = MigBandsInput.AddNewMigBandAfterError(currentPopArray, ancestralPopArray, migBandNewInput, 1, textFields);
//            } else if (!nameOfPopExists(textFields[0].getText(), currentPopArray, ancestralPopArray)) {
//                JOptionPane.showMessageDialog(null, "No population with the given source exists", "message", JOptionPane.ERROR_MESSAGE);                  
//                onComplete = MigBandsInput.AddNewMigBandAfterError(currentPopArray, ancestralPopArray, migBandNewInput, 1, textFields);
//            } else if (!nameOfPopExists(textFields[1].getText(), currentPopArray, ancestralPopArray)) {
//                JOptionPane.showMessageDialog(null, "No population with the given target exists", "message", JOptionPane.ERROR_MESSAGE);
//                onComplete = MigBandsInput.AddNewMigBandAfterError(currentPopArray, ancestralPopArray, migBandNewInput, 1, textFields);
//            } else if ((textFields[0].getText().equals(textFields[1].getText()))) {
//                JOptionPane.showMessageDialog(null, "source and target can't be the same population", "message", JOptionPane.ERROR_MESSAGE);
//                onComplete = MigBandsInput.AddNewMigBandAfterError(currentPopArray, ancestralPopArray, migBandNewInput, 1, textFields);
//            } else if (sourceTargetPairingExists(textFields[0].getText(), textFields[1].getText())){
//                JOptionPane.showMessageDialog(null, "A migration band from " + textFields[0].getText() + " to " + textFields[1].getText() + " already exists" , "message", JOptionPane.ERROR_MESSAGE);
//                onComplete = MigBandsInput.AddNewMigBandAfterError(currentPopArray, ancestralPopArray, migBandNewInput, 1, textFields);
//            } else if (!popIsDaughterOfAncestral(textFields[0].getText(), textFields[1].getText(), currentPopArray, ancestralPopArray)) {
//                JOptionPane.showMessageDialog(null, "the chosen populations can't be a daughter and an ancestral", "message", JOptionPane.ERROR_MESSAGE); 
//                onComplete = MigBandsInput.AddNewMigBandAfterError(currentPopArray, ancestralPopArray, migBandNewInput, 1, textFields);
//            } else {    
//                for (int i = 0; i < NUM_OF_FIELDS; i++) {
//                    migBandInput[i] = textFields[i].getText();
//                }
//                migBandNewInput[0] = migBandInput[0];
//                migBandNewInput[1] = migBandInput[1];
//                onComplete = 2;
//            }
//        }
//        else {
//            onComplete = 0;
//            migBandNewInput[0] = "";
//            migBandNewInput[1] = "";
//        }
//        return onComplete;
//    }
    
    private static int AddNewMigBandAfterError(BSTNode[] currentPopArray, BSTNode[] ancestralPopArray, String[] migBandNewInput, int onComplete, JTextField[] textFields)
    {
     JLabel label;
        String[] migBandInput = {"", ""};
        JPanel addMigBandPanel = new JPanel(new SpringLayout());
        
        
        label = new JLabel(SOURCE, JLabel.TRAILING);
        addMigBandPanel.add(label);
        label.setLabelFor(textFields[0]);
        addMigBandPanel.add(textFields[0]);
        label = new JLabel(TARGET, JLabel.TRAILING);
        addMigBandPanel.add(label);
        label.setLabelFor(textFields[1]);
        addMigBandPanel.add(textFields[1]);             
        

        SpringUtilities.makeCompactGrid(addMigBandPanel,
                
                
                1, NUM_OF_FIELDS * 2, //rows, cols
                6, 6, //initX, initY
                6, 6);       //xPad, yPad

        int result = JOptionPane.showConfirmDialog(null, addMigBandPanel, "Insert the source, target and mig-rate-print of the Migration-Band", JOptionPane.OK_CANCEL_OPTION);

        if (result == JOptionPane.OK_OPTION) {
            if ("".equals(textFields[0].getText()) || "".equals(textFields[1].getText())) {
                JOptionPane.showMessageDialog(null, "One or more of the fields is empty", "message", JOptionPane.ERROR_MESSAGE);
                onComplete = MigBandsInput.AddNewMigBandAfterError(currentPopArray, ancestralPopArray, migBandNewInput, 1, textFields);
            } else if (!nameOfPopExists(textFields[0].getText(), currentPopArray, ancestralPopArray)) {
                JOptionPane.showMessageDialog(null, "No population with the given source exists", "message", JOptionPane.ERROR_MESSAGE);                  
                onComplete = MigBandsInput.AddNewMigBandAfterError(currentPopArray, ancestralPopArray, migBandNewInput, 1, textFields);
            } else if (!nameOfPopExists(textFields[1].getText(), currentPopArray, ancestralPopArray)) {
                JOptionPane.showMessageDialog(null, "No population with the given target exists", "message", JOptionPane.ERROR_MESSAGE);
                onComplete = MigBandsInput.AddNewMigBandAfterError(currentPopArray, ancestralPopArray, migBandNewInput, 1, textFields);
            } else if ((textFields[0].getText().equals(textFields[1].getText()))) {
                JOptionPane.showMessageDialog(null, "source and target can't be the same population", "message", JOptionPane.ERROR_MESSAGE);
                onComplete = MigBandsInput.AddNewMigBandAfterError(currentPopArray, ancestralPopArray, migBandNewInput, 1, textFields);
            } else if (sourceTargetPairingExists(textFields[0].getText(), textFields[1].getText())){
                JOptionPane.showMessageDialog(null, "A migration band from " + textFields[0].getText() + " to " + textFields[1].getText() + " already exists" , "message", JOptionPane.ERROR_MESSAGE);
                onComplete = MigBandsInput.AddNewMigBandAfterError(currentPopArray, ancestralPopArray, migBandNewInput, 1, textFields);
            } else if (!popIsDaughterOfAncestral(textFields[0].getText(), textFields[1].getText(), currentPopArray, ancestralPopArray)) {
                JOptionPane.showMessageDialog(null, "the chosen populations can't be a daughter and an ancestral", "message", JOptionPane.ERROR_MESSAGE); 
                onComplete = MigBandsInput.AddNewMigBandAfterError(currentPopArray, ancestralPopArray, migBandNewInput, 1, textFields);
            } else {    
                for (int i = 0; i < NUM_OF_FIELDS; i++) {
                    migBandInput[i] = textFields[i].getText();
                }
                migBandNewInput[0] = migBandInput[0];
                migBandNewInput[1] = migBandInput[1];
                onComplete = 2;
            }
        }
        else {
            onComplete = 0;
            migBandNewInput[0] = "";
            migBandNewInput[1] = "";
        }
        return onComplete;   
    }
    
    
    private static boolean nameOfPopExists(String name, BSTNode[] currentPopArray, BSTNode[] ancestralPopArray) {
        for(int i = 0; i < currentPopArray.length; i++) {
            if (name.equals(currentPopArray[i].data))
                return true;
        }
        for(int i = 0; i < ancestralPopArray.length; i++) {
            if (name.equals(ancestralPopArray[i].data))
                return true;
        }
        
        return false;
    }
    
    private static boolean sourceTargetPairingExists(String source, String target) {
        int size = ControlFileGeneratorUI.migBandsSourceList.size();
        for(int i = 0; i < size; i++) {
            if(ControlFileGeneratorUI.migBandsSourceList.get(i).equals(source)) {
                if(ControlFileGeneratorUI.migBandsTargetList.get(i).equals(target))
                    return true;
            }
        }
        return false;
    }
    
    private static boolean popIsDaughterOfAncestral(String source, String target, BSTNode[] currentPopArray, BSTNode[] ancestralPopArray) {
        int indexCurSource = -1;
        int indexCurTarget = -1;
        int indexAncSource = -1;
        int indexAncTarget = -1;
        for(int i = 0; i < currentPopArray.length; i++) {
            if(source.equals(currentPopArray[i].data))
                indexCurSource = i;
            if(target.equals(currentPopArray[i].data))
                indexCurTarget = i;
        }
        for(int i = 0; i < ancestralPopArray.length; i++) {
            if(source.equals(ancestralPopArray[i].data))
                indexAncSource = i;
            if(target.equals(ancestralPopArray[i].data))
                indexAncTarget = i;
        }
        if(indexCurSource != -1 && indexCurTarget != -1)
            //both are current pop
            return true;
        else if(indexCurSource != -1 && indexAncTarget != -1)
            return checkRelations(indexCurSource, indexAncTarget, currentPopArray, ancestralPopArray);
        else if(indexAncSource != -1 && indexCurTarget != -1)
            return checkRelations(indexAncSource, indexCurTarget, ancestralPopArray, currentPopArray);
        else if(indexAncSource != -1 && indexAncTarget != -1)
            return checkRelations(indexAncSource, indexAncTarget, ancestralPopArray, ancestralPopArray);
        else
            return false;
    }
    
    private static boolean checkRelations(int indexSource, int indexTarget, BSTNode[] source, BSTNode[] target) {
        String sourceData = source[indexSource].data;
        String targetData = target[indexTarget].data;
        BSTNode sourceTraversal = source[indexSource];
        BSTNode targetTraversal = target[indexTarget];
        while(sourceTraversal != null) {
            if(targetData.equals(sourceTraversal.data))
                return false;
            sourceTraversal = sourceTraversal.parent;
        }
        while(targetTraversal != null) {
            if(sourceData.equals(targetTraversal.data))
                return false;
            targetTraversal = targetTraversal.parent;
        }
        return true;
    }
    
    public static void DeleteMigBands(int migBandsCounter) {
        JLabel label = new JLabel("", JLabel.TRAILING);
        JTextField[] sourceTextFields = new JTextField[migBandsCounter];
        JTextField[] targetTextFields = new JTextField[migBandsCounter];
        JCheckBox[] markToDeleteButtons = new JCheckBox[migBandsCounter];
        String[] newMigBandsSource = new String[migBandsCounter];
        String[] newMigBandsTarget = new String[migBandsCounter];

        for (int i = 0; i < migBandsCounter; i++) {
            sourceTextFields[i] = new JTextField(5);
            sourceTextFields[i].setEnabled(false);
            sourceTextFields[i].setText(ControlFileGeneratorUI.migBandsSourceList.get(i));
            sourceTextFields[i].setDisabledTextColor(Color.black);
            targetTextFields[i] = new JTextField(5);
            targetTextFields[i].setEnabled(false);
            targetTextFields[i].setText(ControlFileGeneratorUI.migBandsTargetList.get(i));
            targetTextFields[i].setDisabledTextColor(Color.black);
            markToDeleteButtons[i] = new JCheckBox();
            markToDeleteButtons[i].setSelected(false);
            
            newMigBandsSource[i] = "";
            newMigBandsTarget[i] = "";
        }

        JPanel modifyMigBandPanel = new JPanel(new SpringLayout());

        for (int i = 0; i < migBandsCounter; i++) {
            label.setLabelFor(markToDeleteButtons[i]);
            modifyMigBandPanel.add(markToDeleteButtons[i]);
            label = new JLabel(SOURCE, JLabel.TRAILING);
            modifyMigBandPanel.add(label);
            label.setLabelFor(sourceTextFields[i]);
            modifyMigBandPanel.add(sourceTextFields[i]);
            label = new JLabel(TARGET, JLabel.TRAILING);
            modifyMigBandPanel.add(label);
            label.setLabelFor(targetTextFields[i]);
            modifyMigBandPanel.add(targetTextFields[i]);
        }

        SpringUtilities.makeCompactGrid(modifyMigBandPanel,
                migBandsCounter, 5, //rows, cols
                6, 6, //initX, initY
                6, 6);       //xPad, yPad

        int result = JOptionPane.showConfirmDialog(null, modifyMigBandPanel, "Mark the Migration-Bands you wish to delete", JOptionPane.OK_CANCEL_OPTION);
        
        if (result == JOptionPane.OK_OPTION) {
            for(int i = 0; i < migBandsCounter; i++) {
                if(!markToDeleteButtons[i].isSelected()) {
                    newMigBandsSource[i] = sourceTextFields[i].getText();
                    newMigBandsTarget[i] = targetTextFields[i].getText();
                }
            }
            ControlFileGeneratorUI.migBandsSourceList.clear();
            ControlFileGeneratorUI.migBandsTargetList.clear();
            
            for(int i = 0; i < migBandsCounter; i++) {
                if(!"".equals(newMigBandsSource[i])) {
                    ControlFileGeneratorUI.migBandsSourceList.add(newMigBandsSource[i]);
                    ControlFileGeneratorUI.migBandsTargetList.add(newMigBandsTarget[i]);
                }
            }
        }
    }
}