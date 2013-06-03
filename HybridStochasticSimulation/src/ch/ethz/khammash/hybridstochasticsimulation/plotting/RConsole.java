package ch.ethz.khammash.hybridstochasticsimulation.plotting;

import java.awt.Frame;
import java.io.BufferedReader;
import java.io.InputStreamReader;

import javax.swing.JFileChooser;
import javax.swing.filechooser.FileNameExtensionFilter;

import org.rosuda.JRI.RMainLoopCallbacks;
import org.rosuda.JRI.Rengine;

class RConsole implements RMainLoopCallbacks {

	Frame parent;

	public RConsole(Frame parent) {
		this.parent = parent;
	}

    public void rWriteConsole(Rengine re, String text, int oType) {
        System.out.print(text);
    }

    public void rBusy(Rengine re, int which) {
        System.out.println("rBusy(" + which + ")");
    }

    public String rReadConsole(Rengine re, String prompt, int addToHistory) {
        System.out.print(prompt);
        try {
            BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
            String s = br.readLine();
            return (s == null || s.length() == 0) ? s : s + "\n";
        } catch (Exception e) {
            System.out.println("jriReadConsole exception: " + e.getMessage());
        }
        return null;
    }

    public void rShowMessage(Rengine re, String message) {
        System.out.println("rShowMessage \"" + message + "\"");
    }

    public String rChooseFile(Rengine re, int newFile) {
		JFileChooser fc = new CustomJFileChooser();
		fc.setFileFilter(new FileNameExtensionFilter("PDF documents", "pdf"));
		int returnValue = fc.showSaveDialog(parent);
		String result = null;
		if (returnValue == JFileChooser.APPROVE_OPTION)
			result = fc.getSelectedFile().getAbsolutePath();
		return result;
    }

    public void   rFlushConsole (Rengine re) {
    }

    public void   rLoadHistory  (Rengine re, String filename) {
    }			

    public void   rSaveHistory  (Rengine re, String filename) {
    }

}
