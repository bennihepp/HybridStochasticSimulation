package ch.ethz.khammash.hybridstochasticsimulation.plotting;

import java.awt.Font;
import java.io.ByteArrayOutputStream;
import java.io.InputStream;

import com.itextpdf.awt.DefaultFontMapper;
import com.itextpdf.text.pdf.BaseFont;

// TODO
public class PdfFontMapper extends DefaultFontMapper {

	@Override
	public BaseFont awtToPdf(Font font) {
		System.out.println("Mapping font: "+font.getSize2D());
		BaseFont bf = super.awtToPdf(font.deriveFont(font.getSize2D() * 500));
		return bf;
	}

//	public BaseFont awtToPdf(Font font) {
//
//		try {
//			if(font.getFamily().equalsIgnoreCase("Verdana")) {
//				if(font.isBold()) {
//					if(font.isItalic()) {
//						return this.getBaseFontFromFile("/META-INF/fonts/verdana/", "VERDANAZ.TTF");
//					}
//					return this.getBaseFontFromFile("/META-INF/fonts/verdana/", "VERDANAB.TTF");
//				} else if(font.isItalic()) {
//					return this.getBaseFontFromFile("/META-INF/fonts/verdana/", "VERDANAI.TTF");
//				} else {
//					return this.getBaseFontFromFile("/META-INF/fonts/verdana/", "VERDANA.TTF");
//				}
//			} else { //Times new Roman is default
//				if(font.isBold()) {
//					if(font.isItalic()) {
//						return this.getBaseFontFromFile("/META-INF/fonts/timesnewroman/", "TIMESBI.TTF");
//					}
//					return this.getBaseFontFromFile("/META-INF/fonts/timesnewroman/", "TIMESBD.TTF");
//				} else if(font.isItalic()) {
//					return this.getBaseFontFromFile("/META-INF/fonts/timesnewroman/", "TIMESI.TTF");
//				} else {
//					return this.getBaseFontFromFile("/META-INF/fonts/timesnewroman/", "TIMES.TTF");
//				}
//			}
//			
//		} catch(Exception e) {
//			e.printStackTrace();
//			throw new RuntimeException(e);
//		}
//	}

	public Font pdfToAwt(BaseFont baseFont, int size) {
		throw new UnsupportedOperationException();
	}

	/**
	 * To get a {@link BaseFont} from a file on the filesystem or in a jar. See: http://www.mail-archive.com/itext-questions@lists.sourceforge.net/msg02691.html
	 *
	 * @param directory
	 * @param filename
	 * @return
	 * @throws Exception
	 */
	private BaseFont getBaseFontFromFile(String directory, String filename) throws Exception {
	
		InputStream is = null;
		try {
			is = PdfFontMapper.class.getResourceAsStream(directory + filename);

			ByteArrayOutputStream bos = new ByteArrayOutputStream();
			byte[] buf = new byte[1024];
			while(true) {
				int size = is.read(buf);
				if(size < 0) {
					break;
				}
				bos.write(buf, 0, size);
			}
			buf = bos.toByteArray();
			BaseFont bf = BaseFont.createFont(filename, BaseFont.WINANSI, BaseFont.NOT_EMBEDDED, BaseFont.NOT_CACHED, buf, null);
			return bf;
		} finally {
			if(is != null) {
				is.close();
			}
		}
	}

}
