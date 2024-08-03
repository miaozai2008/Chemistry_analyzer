#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QDebug>
#include <QLabel>
#include <QMessageBox>
#include <QTime>

MainWindow::MainWindow(QWidget* parent)
	: QMainWindow(parent)
	, ui(new Ui::MainWindow)
{
	//配置界面
	ui->setupUi(this);
	ui->statusBar->showMessage(QStringLiteral("欢迎使用热化学反应解析器"), 8000);
	ui->statusBar->addPermanentWidget(new QLabel("DYZ(C) v5.0 2024.7.13", this));
	//连接core
	core->setProgram("Core.exe");
	connect(core, &QProcess::errorOccurred, this, [=]() {
		ui->textBrowser->setHtml(QStringLiteral("|Core进程错误:"));
		ui->textBrowser->append(core->errorString());
		ui->textBrowser->append(QStringLiteral("请在Core.exe就位后重启软件"));
		ui->pushButton->setDisabled(true);
		});
	connect(core, &QProcess::readyReadStandardOutput, this, [=]() {
		static QStringDecoder toUtf = QStringDecoder(QStringDecoder::System);
		ui->textBrowser->append(toUtf(core->readAllStandardOutput()));
		});
	core->start(QIODevice::ReadWrite);
}

MainWindow::~MainWindow() { core->write("[exit]\n"); core->disconnect(); delete ui; }

void MainWindow::on_action_this_triggered() {//解析器介绍
	ui->textBrowser->setHtml(QStringLiteral("|关于Chemistry Analyzer:"));
	ui->textBrowser->append(QStringLiteral("版本:v5.0 2024.7.13\n\n(人名皆为首字母)"));
	ui->textBrowser->append(QStringLiteral("指导老师:ZX\n版权:DYZ(C)"));
	ui->textBrowser->append(QStringLiteral("开发工具:VS2022 Enterprise, Qt6.7.2"));
	ui->textBrowser->append(QStringLiteral("团队成员:(字母升序)\nCYX DYZ GJC PJ WYS WZX YHX YZT ZHY"));
}

void MainWindow::on_action_H_triggered() {//使用辅助
	ui->textBrowser->setHtml(QStringLiteral("|使用辅助"));
	ui->textBrowser->append(QStringLiteral("合法物质:e<-> NH4<+>(aq) 1/2O2 [KAl(SO4)2·12H2O] 支持括号嵌套"));
	ui->textBrowser->append(QStringLiteral("输入示范:1/3H2O(l)+Cl2(g)--H<+>(aq)+Cl<->(aq)+HClO(aq) 空格不限"));
	ui->textBrowser->append(QStringLiteral("STP(标准温度&气压):25°C,101kPa"));
}

void MainWindow::on_action_Qt_triggered() { QMessageBox::aboutQt(this, "About Qt"); }//关于Qt

void MainWindow::on_action_code_triggered() {//github源代码
	ui->textBrowser->setHtml(QStringLiteral("|github项目地址:<a href="
		"\"https://github.com/miaozai2008/Chemistry-analyzer\">右键点击前往</a>"));
}

void MainWindow::on_action_issue_triggered() {//提交issue
	ui->textBrowser->setHtml(QStringLiteral("|github项目issue:<a href="
		"\"https://github.com/miaozai2008/Chemistry-analyzer/issues\">右键点击前往</a>"));
}

void MainWindow::on_pushButton_clicked() {//点击按钮
	ui->textBrowser->setHtml(QTime::currentTime().toString("|hh:mm:ss"));
	if (ui->lineEdit->text().isEmpty() || ui->lineEdit->text() == "[exit]") {
		ui->textBrowser->append(QStringLiteral("化学式错误"));
		return;
	}
	core->write((ui->lineEdit->text().trimmed().remove(" ") + "\n").toLocal8Bit());
	core->write((QString::number(ui->doubleSpinBoxT->value() + 273.15) + "\n").toLocal8Bit());
	core->write((QString::number(ui->doubleSpinBoxP->value()) + "\n").toLocal8Bit());
}

void MainWindow::error() {
	ui->textBrowser->setHtml(QStringLiteral("|Core进程错误:"));
	ui->textBrowser->append(core->errorString());
	ui->textBrowser->append(QStringLiteral("请在Core.exe就位后重启软件"));
	QMessageBox::critical(this, QStringLiteral("进程错误"), QStringLiteral("Core进程出现错误"));
	ui->pushButton->setDisabled(true);
}

void MainWindow::output() {
	QByteArray arr = core->readAllStandardOutput();
	static QStringDecoder toUtf = QStringDecoder(QStringDecoder::System);
	ui->textBrowser->append(toUtf(arr));
}