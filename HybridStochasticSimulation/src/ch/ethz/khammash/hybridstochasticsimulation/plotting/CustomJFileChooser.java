package ch.ethz.khammash.hybridstochasticsimulation.plotting;

import java.io.File;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

public class CustomJFileChooser extends JFileChooser {

	private static final long serialVersionUID = -1906630961401909120L;

	public void approveSelection() {
        File f = getSelectedFile();
        if (f.exists() && getDialogType() == SAVE_DIALOG) {
			int returnValue = JOptionPane.showConfirmDialog(this,
					"The file exists, overwrite?", "Existing file",
					JOptionPane.YES_NO_CANCEL_OPTION);
            switch (returnValue) {
            case JOptionPane.YES_OPTION:
                super.approveSelection();
                break;
            case JOptionPane.CANCEL_OPTION:
                cancelSelection();
                break;
            default:
            }
        } else
            super.approveSelection();
    }

}
